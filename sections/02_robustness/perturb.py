# Builds the perturbed reference for Results 2.2 (robustness).
#
# Pair 0: scRNA 01_mouse-ssp (5665 x 45768) x ST 01_mouse-ssp (osmFISH, 4839 x 33).
# The perturbation touches the reference only; the slice is untouched, so the
# correct behaviour is known a priori.
#
#   split.h5ad      obs["cellType_split"]  4 types
#       Glutamatergic is halved into _A / _B by pairing each cell with its
#       nearest neighbour in panel space and sending one of each pair to each
#       half, so the two halves are alike in everything the slice measures. The
#       ten flattest non-panel genes inside Glutamatergic then receive, in the B
#       half only, the count distribution of one of the ten strongest
#       Glutamatergic markers (resampled, jittered, shuffled). The 33 shared
#       genes are never touched, so the split leaves no trace in G n M.
#       Everything here uses the coarse cellType annotation only.
#
# Usage:
#   & ".../spatialaim_env/python.exe" perturb.py

import json
import logging

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

from spatialaim_experiments import paths

logger = logging.getLogger(__name__)

SC_PATH = paths.DATASETS / "scRNA" / "01_mouse-ssp.h5ad"
ST_PATH = paths.DATASETS / "ST" / "01_mouse-ssp.h5ad"
# the perturbed reference is a dataset, so it is written where the datasets live
OUT_DIR = paths.DATASETS / "Robustness"

CELL_TYPE_KEY = "cellType"
SPLIT_KEY = "cellType_split"

SPLIT_TARGET = "Glutamatergic"
SPLIT_SEED = 0
N_MARKER_GENES = 10
NOISE_SIGMA = 0.3  # multiplicative lognormal jitter on the resampled counts
# Donor markers: how strongly a gene must be on inside the type it marks.
MIN_MARKER_DETECTION = 0.5
# How many neighbours to look at when matching the target cells into pairs.
PAIRING_NEIGHBOURS = 30


def dense(x) -> np.ndarray:
    return np.asarray(x.todense()) if sp.issparse(x) else np.asarray(x)


def matched_halves(
    panel_counts: np.ndarray,
    target: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Halve ``target`` into two panel-matched halves.

    Each target cell is paired with its nearest unpaired neighbour in raw panel
    space, and one of each pair goes to A, the other to B. A plain coin flip
    would leave the two halves differing by sampling noise, and in this panel
    that noise is largest exactly in the layer markers -- enough for a mapper to
    sort the spots by cortical layer. Pairing replaces the within-type spread by
    the within-pair spread, which is far smaller, so the two halves look alike to
    anything that only sees the panel.
    """
    from sklearn.neighbors import NearestNeighbors

    G = panel_counts[target]
    n = len(G)
    k = min(PAIRING_NEIGHBOURS, n - 1)
    nbr = NearestNeighbors(n_neighbors=k + 1).fit(G).kneighbors(G, return_distance=False)

    taken = np.zeros(n, dtype=bool)
    pairs: list[tuple[int, int]] = []
    for i in rng.permutation(n):
        if taken[i]:
            continue
        for j in nbr[i][1:]:
            if not taken[j]:
                taken[i] = taken[j] = True
                pairs.append((int(i), int(j)))
                break
    # cells whose whole neighbourhood was already taken, plus the odd one out
    left = np.where(~taken)[0]
    rng.shuffle(left)
    pairs += [(int(left[t]), int(left[t + 1])) for t in range(0, len(left) - 1, 2)]
    pair_idx = np.array(pairs)

    flip = rng.random(len(pair_idx)) < 0.5
    a = np.where(flip, pair_idx[:, 0], pair_idx[:, 1])
    b = np.where(flip, pair_idx[:, 1], pair_idx[:, 0])
    if len(left) % 2:  # one cell cannot be paired; it joins a half at random
        odd = int(left[-1])
        if rng.random() < 0.5:
            a = np.append(a, odd)
        else:
            b = np.append(b, odd)

    within = np.linalg.norm(G[pair_idx[:, 0]] - G[pair_idx[:, 1]], axis=1)
    shuffled = rng.permutation(n)
    across = np.linalg.norm(G[shuffled[: len(pair_idx)]] - G[shuffled[-len(pair_idx):]], axis=1)
    ca, cb = G[a].mean(0), G[b].mean(0)
    sem = np.sqrt(2) * G.std(0) / np.sqrt(n // 2)  # sd of a coin-flip half-mean difference
    stats = {
        "n_pairs": int(len(pair_idx)),
        "n_unpairable": int(len(left)),
        "median_within_pair_distance": float(np.median(within)),
        "median_random_pair_distance": float(np.median(across)),
        "panel_centroid_cosine": float(ca @ cb / (np.linalg.norm(ca) * np.linalg.norm(cb))),
        "max_gene_imbalance_in_coinflip_sd": float(
            np.max(np.abs(ca - cb) / np.where(sem > 0, sem, np.inf))
        ),
    }
    return np.sort(target[a]), np.sort(target[b]), stats


def pick_recipient_genes(
    X: np.ndarray,
    var_names: pd.Index,
    panel: set[str],
    target: np.ndarray,
    rng: np.random.Generator,
    n: int,
) -> list[str]:
    """The n non-panel genes with the lowest count variance inside the split type.

    Lying outside the panel is what makes the injected signature invisible to the
    slice; the lowest variance inside the target cells is what makes overwriting
    half of them destroy no signal. The lowest variance is heavily tied (a gene
    silent in every cell of the type has variance 0), so genes below the cut-off
    are taken outright and the remaining slots are drawn at random from the ones
    sitting exactly on it, rather than letting the alphabet decide.
    """
    cand = np.where(np.array([g not in panel for g in var_names]))[0]
    v = X[np.ix_(target, cand)].var(axis=0)
    cut = np.sort(v)[n - 1]
    below, tied = cand[v < cut], cand[v == cut]
    logger.info(
        "  %d non-panel genes, %d below the cut-off variance %.4g, %d tied on it",
        len(cand), len(below), cut, len(tied),
    )
    drawn = rng.choice(tied, size=n - len(below), replace=False)
    return [str(var_names[j]) for j in np.sort(np.concatenate([below, drawn]))]


def pick_donor_markers(
    X: np.ndarray,
    var_names: pd.Index,
    cell_type: np.ndarray,
    donor_type: str,
    exclude: set[str],
    n: int,
) -> list[str]:
    """The n strongest markers of ``donor_type`` under the coarse annotation.

    Ranked by the Wilcoxon effect size (AUC) of the type against all other
    cells, requiring the gene to be on in at least ``MIN_MARKER_DETECTION`` of
    the type's own cells. The U statistic is read off one global rank transform
    (rank sums per group) rather than a test per gene.
    """
    usable = np.where(np.array([g not in exclude for g in var_names]))[0]
    idx = np.where(cell_type == donor_type)[0]
    m, r = len(idx), X.shape[0] - len(idx)

    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    ranks = rankdata(np.log1p(X[:, usable] / lib * 1e4), axis=0)
    auc = (ranks[idx].sum(axis=0) - m * (m + 1) / 2.0) / (m * r)

    on = (X[np.ix_(idx, usable)] > 0).mean(axis=0)
    score = np.where(on >= MIN_MARKER_DETECTION, auc, -np.inf)
    logger.info("  %d candidate marker genes, %d on in >= %.0f%% of %s cells",
                len(usable), int(np.isfinite(score).sum()), 100 * MIN_MARKER_DETECTION, donor_type)
    return [str(var_names[usable[j]]) for j in np.argsort(score)[::-1][:n]]


def build_split(sc: ad.AnnData, panel: set[str]) -> tuple[ad.AnnData, dict]:
    X = dense(sc.X).astype(np.float32)
    # separate streams so that changing the halving does not reshuffle the genes
    rng_split = np.random.default_rng(SPLIT_SEED)
    rng_genes = np.random.default_rng(SPLIT_SEED + 1000)
    rng_write = np.random.default_rng(SPLIT_SEED + 2000)

    labels = sc.obs[CELL_TYPE_KEY].astype(str).to_numpy()
    target = np.where(labels == SPLIT_TARGET)[0]
    panel_cols = [g for g in sc.var_names if g in panel]
    panel_counts = np.asarray(sc[:, panel_cols].X, dtype=float)

    logger.info("Pairing %s cells in panel space ...", SPLIT_TARGET)
    idx_a, idx_b, pairing = matched_halves(panel_counts, target, rng_split)
    logger.info("Splitting %s: %d -> A %d / B %d", SPLIT_TARGET, len(target), len(idx_a), len(idx_b))
    logger.info("  %d pairs, %d unpairable; within-pair distance %.1f vs %.1f for random pairs",
                pairing["n_pairs"], pairing["n_unpairable"],
                pairing["median_within_pair_distance"], pairing["median_random_pair_distance"])
    logger.info("  panel centroid cosine A vs B %.7f; largest per-gene imbalance %.2f coin-flip sd",
                pairing["panel_centroid_cosine"], pairing["max_gene_imbalance_in_coinflip_sd"])

    logger.info("Choosing recipient genes ...")
    recipients = pick_recipient_genes(X, sc.var_names, panel, target, rng_genes, N_MARKER_GENES)
    logger.info("  recipients: %s", recipients)

    logger.info("Choosing donor markers of %s ...", SPLIT_TARGET)
    donors = pick_donor_markers(X, sc.var_names, labels, SPLIT_TARGET, set(recipients) | panel, N_MARKER_GENES)
    logger.info("  donors: %s", donors)

    pairs = []
    for donor_gene, recipient in zip(donors, recipients):
        d_col = sc.var_names.get_loc(donor_gene)
        r_col = sc.var_names.get_loc(recipient)
        pool = X[target, d_col]  # the donor's counts inside the type being split
        drawn = rng_write.choice(pool, size=len(idx_b), replace=True)
        drawn = np.rint(drawn * rng_write.lognormal(0.0, NOISE_SIGMA, size=len(idx_b)))
        before = float(X[idx_b, r_col].mean())
        X[idx_b, r_col] = np.maximum(drawn, 0.0).astype(np.float32)
        pairs.append(
            {
                "recipient": recipient,
                "recipient_mean_counts_before": before,
                "donor": donor_gene,
                "donor_mean_counts": float(pool.mean()),
                "written_mean_counts": float(drawn.mean()),
            }
        )

    new_labels = labels.copy()
    new_labels[idx_a] = f"{SPLIT_TARGET}_A"
    new_labels[idx_b] = f"{SPLIT_TARGET}_B"

    obs = sc.obs.copy()
    obs[SPLIT_KEY] = pd.Categorical(new_labels)
    half_col = np.array([""] * sc.n_obs, dtype=object)
    half_col[idx_a] = "A"
    half_col[idx_b] = "B"
    obs["split_half"] = pd.Categorical(half_col)

    out = ad.AnnData(X=X, obs=obs, var=sc.var.copy())
    meta = {
        "type": "invisible_split",
        "target_type": SPLIT_TARGET,
        "seed": SPLIT_SEED,
        "n_target": int(len(target)),
        "n_a": int(len(idx_a)),
        "n_b": int(len(idx_b)),
        "n_marker_genes": N_MARKER_GENES,
        "noise_sigma_lognormal": NOISE_SIGMA,
        "halving_rule": (
            f"nearest-neighbour pairing in raw panel space ({PAIRING_NEIGHBOURS} candidate "
            f"neighbours), one of each pair to each half, orientation by coin flip"
        ),
        "pairing": pairing,
        "recipient_rule": (
            f"lowest count variance within {SPLIT_TARGET}, among non-panel genes; "
            f"ties on the cut-off broken by a random draw, seed {SPLIT_SEED}"
        ),
        "donor_rule": (
            f"strongest {SPLIT_TARGET} markers by Wilcoxon AUC vs all other cells, "
            f"on in >= {MIN_MARKER_DETECTION:.0%} of {SPLIT_TARGET} cells; counts "
            f"resampled from within {SPLIT_TARGET}"
        ),
        "annotation_used": CELL_TYPE_KEY,
        "gene_pairs": pairs,
        "annotation_column": SPLIT_KEY,
    }
    # h5ad cannot store nested/mixed structures in uns -> keep it as JSON text;
    # manifest.json carries the structured version.
    out.uns["perturbation"] = json.dumps(meta)
    return out, meta


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = OUT_DIR / "manifest.json"

    st = ad.read_h5ad(ST_PATH, backed="r")
    panel = set(st.var_names)
    logger.info("Panel: %d genes", len(panel))

    logger.info("Loading %s ...", SC_PATH)
    sc = ad.read_h5ad(SC_PATH)
    logger.info("  %d cells x %d genes", sc.n_obs, sc.n_vars)

    manifest: dict = {}
    manifest["pair"] = {"sc": SC_PATH.name, "st": ST_PATH.name, "shared_genes": len(set(sc.var_names) & panel)}
    manifest["clean"] = {"path": str(SC_PATH), "annotation_column": CELL_TYPE_KEY}

    split, meta = build_split(sc, panel)
    path = OUT_DIR / "split.h5ad"
    split.write_h5ad(path)
    meta["path"] = str(path)
    manifest["split"] = meta
    logger.info("Wrote %s", path)

    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    logger.info("Wrote %s", manifest_path)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
