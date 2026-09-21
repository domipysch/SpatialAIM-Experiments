# Numbers for the invisible-split experiment (Results 2.2, first half).
#
# Three questions, in the order the subsection argues them:
#   1. Is the split really invisible in the panel, while being a real cell type in
#      the reference?
#   2. Handed the perturbed annotation, what does the mapper return at k_004 --
#      how much of the tissue goes to each phantom half, and does the division
#      carry any spatial structure?
#   3. What does SpatialAIM's level selection do with the same input, and does the level
#      it proposes reproduce the unperturbed mapping?
#
# One mapper per invocation (``--mapper``, default Tangram). Writes
# Data/02_Results/robustness/invisible_split/_analysis/<mapper>/
#   results.json   everything
#   metrics.csv    the handful of numbers the text quotes
#
# Run from spatialaim_env (see the section README) -- imports squidpy via spatialaim.

import argparse
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import adjusted_rand_score

from spatialaim.analysis.topology import spatial_connectivities
from spatialaim.metrics.topology import local_spatial_purity, permutation_test

from spatialaim_experiments import paths
from spatialaim_experiments.style import MAPPERS

logger = logging.getLogger(__name__)

ROOT = paths.ROBUSTNESS / "invisible_split"
ST_PATH = paths.DATASETS / "ST" / "01_mouse-ssp.h5ad"
SC_ORIGINAL = paths.DATASETS / "scRNA" / "01_mouse-ssp.h5ad"
SC_SPLIT = paths.DATASETS / "Robustness" / "split.h5ad"

MAPPER = "tangram"
TARGET = "Glutamatergic"
A, B = f"{TARGET}_A", f"{TARGET}_B"
K_IDENTITY, K_BELOW = 4, 3
N_PERM = 200


class Run:
    """One sweep on disk (``<reference>/<mapper>/``)."""

    def __init__(self, reference: str):
        self.reference = reference
        self.dir = ROOT / reference / MAPPER
        self.table = pd.read_csv(self.dir / "k_comparison.csv")
        self.selection = json.loads((self.dir / "k_selection.json").read_text(encoding="utf-8"))
        sc = ad.read_h5ad(self.dir / "start_clustering.h5ad")
        self.start_cluster = sc.obs["start_cluster"].to_numpy().astype(int)
        self.start_cluster_name = sc.obs["start_cluster_name"].astype(str).to_numpy()

    def cut(self, k: int) -> np.ndarray:
        df = pd.read_csv(self.dir / f"k_{k:03d}" / "start_cluster_to_state.csv")
        out = np.zeros(int(self.start_cluster.max()) + 1, dtype=int)
        out[df["start_cluster"].to_numpy()] = df["state"].to_numpy()
        return out

    def P(self, k: int) -> np.ndarray:
        path = self.dir / f"k_{k:03d}" / "spot_to_state_mapping_soft.h5ad"
        return np.asarray(ad.read_h5ad(path).X, dtype=float)

    def hard(self, k: int) -> np.ndarray:
        return np.argmax(self.P(k), axis=1)

    def state_names(self, k: int, undo_split: bool = False) -> list[str]:
        """The annotated types merged into each state, as one readable label.

        With ``undo_split`` the two phantom halves are renamed back to the type
        they came from before the labels are joined, so a state merging exactly
        those two reads as the original type and is comparable with the
        unperturbed run.
        """
        cs = self.cut(k)
        out = []
        for s in range(int(cs.max()) + 1):
            members = set(self.start_cluster_name[np.isin(self.start_cluster, np.where(cs == s)[0])])
            if undo_split:
                members = {TARGET if m in (A, B) else m for m in members}
            out.append(" + ".join(sorted(members)))
        return out

    def state_of(self, k: int, type_name: str) -> int:
        return int(self.cut(k)[self.start_cluster[self.start_cluster_name == type_name][0]])

    def spot_types(self, k: int, undo_split: bool = False) -> np.ndarray:
        return np.array(self.state_names(k, undo_split), dtype=object)[self.hard(k)]

    def score_row(self, k: int) -> dict:
        row = self.table[self.table["k"] == k].iloc[0]
        return {c: float(row[c]) for c in ("recon", "spatial", "coherence", "overall")}

    def best(self) -> dict[str, int | None]:
        return {
            name: (None if hit is None else int(hit["k"]))
            for name, hit in self.selection["best"].items()
        }


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0.0, 1.0)


def invisibility() -> dict:
    """How far apart the two halves are in the panel, and in the whole transcriptome."""
    sc = ad.read_h5ad(SC_SPLIT)
    st = ad.read_h5ad(ST_PATH, backed="r")
    panel = [g for g in sc.var_names if g in set(st.var_names)]
    half = sc.obs["split_half"].astype(str).to_numpy()
    ia, ib = half == "A", half == "B"
    out = {"n_a": int(ia.sum()), "n_b": int(ib.sum())}
    for name, genes in (("panel", panel), ("all_genes", list(sc.var_names))):
        X = np.asarray(sc[:, genes].X, dtype=float)
        ca, cb = X[ia].mean(0), X[ib].mean(0)
        _, p = mannwhitneyu(X[ia], X[ib], axis=0)
        out[name] = {
            "n_genes": len(genes),
            "centroid_cosine": float(ca @ cb / (np.linalg.norm(ca) * np.linalg.norm(cb))),
            "n_significant_q005": int((bh(np.nan_to_num(p, nan=1.0)) < 0.05).sum()),
        }
        del X
    return out


def phantom_division(split: Run, conn) -> dict:
    """What the mapper does with the two phantom halves at the identity level."""
    names = split.state_names(K_IDENTITY)
    sa, sb = split.state_of(K_IDENTITY, A), split.state_of(K_IDENTITY, B)
    P = split.P(K_IDENTITY)
    hard = np.argmax(P, axis=1)
    n = P.shape[0]

    both = np.where((hard == sa) | (hard == sb))[0]
    binary = (hard == sb).astype(int)
    sub = conn[both][:, both]
    observed = local_spatial_purity(binary[both], sub)
    null = permutation_test(
        observed, binary[both], lambda x: local_spatial_purity(x, sub), N_PERM,
        np.random.default_rng(0),
    )
    return {
        "state_names": names,
        "spot_share": {names[s]: float((hard == s).mean()) for s in range(len(names))},
        # Tangram's cluster mode normalises to equal mass per cell type, so only the
        # argmax shares above are comparable; kept here only to make that visible.
        "mass": {names[s]: float(P[:, s].sum() / P.sum()) for s in range(len(names))},
        "n_spots_a": int((hard == sa).sum()),
        "n_spots_b": int((hard == sb).sum()),
        "imbalance": float(abs((hard == sa).sum() - (hard == sb).sum()) / max(len(both), 1)),
        "n_spots_total": int(n),
        "ab_purity": {k: null[k] for k in ("observed", "null_mean", "null_std", "z_score", "p_value")},
    }


def recovery(split: Run, original: Run) -> dict:
    """Does the level below the perturbed annotation reproduce the unperturbed map?"""
    base = original.spot_types(K_BELOW)
    below = split.spot_types(K_BELOW, undo_split=True)
    identity = split.spot_types(K_IDENTITY, undo_split=True)
    return {
        "first_merge_is_phantom_pair": bool(
            split.state_of(K_BELOW, A) == split.state_of(K_BELOW, B)
        ),
        "label_agreement_below": float((below == base).mean()),
        "ari_below": float(adjusted_rand_score(base, below)),
        "label_agreement_identity": float((identity == base).mean()),
        "ari_identity": float(adjusted_rand_score(base, split.spot_types(K_IDENTITY))),
    }


def merge_heights() -> dict:
    """The agglomeration over the four annotated types, on shared-gene centroids."""
    from spatialaim.tree import build_agglomeration_tree

    scaffold = ad.read_h5ad(ROOT / "split" / "reference_scaffold.h5ad", backed="r")
    names = [str(n) for n in scaffold.uns["start_cluster_names"]]
    Z = build_agglomeration_tree(np.asarray(scaffold.uns["start_cluster_centroids_shared"]),
                                 method="average")
    joined, out = [[n] for n in names], []
    for left, right, height, _ in Z:
        members = joined[int(left)] + joined[int(right)]
        joined.append(members)
        out.append({"merge": " + ".join(sorted(members)), "height": float(height)})
    return {"merges": out, "first_merge_height": out[0]["height"],
            "ratio_first_to_second": out[0]["height"] / out[1]["height"]}


def main() -> None:
    global MAPPER
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", default=MAPPER,
                    choices=MAPPERS,
                    help="which mapper's runs to read (default tangram)")
    MAPPER = ap.parse_args().mapper

    out = ROOT / "_analysis" / MAPPER
    out.mkdir(parents=True, exist_ok=True)
    st = ad.read_h5ad(ST_PATH)
    conn = spatial_connectivities(st, 6)

    original, split = Run("original"), Run("split")
    results = {
        "mapper": MAPPER,
        "invisibility": invisibility(),
        "tree": merge_heights(),
        "phantom": phantom_division(split, conn),
        "recovery": recovery(split, original),
        "selection": {
            "original": {"best": original.best(), "pareto_ks": original.selection["pareto_ks"]},
            "split": {"best": split.best(), "pareto_ks": split.selection["pareto_ks"]},
        },
        "scores": {
            "original": original.table.to_dict(orient="list"),
            "split": split.table.to_dict(orient="list"),
        },
        "score_delta_identity_vs_below": {
            key: split.score_row(K_IDENTITY)[key] - split.score_row(K_BELOW)[key]
            for key in ("recon", "spatial", "coherence", "overall")
        },
    }

    inv, ph, rec, sel = (results[k] for k in ("invisibility", "phantom", "recovery", "selection"))
    rows = [
        ("panel genes differing between A and B", f"{inv['panel']['n_significant_q005']} of {inv['panel']['n_genes']}"),
        ("all genes differing between A and B", f"{inv['all_genes']['n_significant_q005']} of {inv['all_genes']['n_genes']}"),
        ("panel centroid cosine A vs B", f"{inv['panel']['centroid_cosine']:.7f}"),
        ("first tree merge (height)", f"{results['tree']['merges'][0]['merge']} at {results['tree']['first_merge_height']:.4f}"),
        ("first merge vs next, as a ratio", f"{results['tree']['ratio_first_to_second']:.4f}"),
        (f"spot share on {A} (k_{K_IDENTITY:03d})", f"{ph['spot_share'][A]:.3f}"),
        (f"spot share on {B} (k_{K_IDENTITY:03d})", f"{ph['spot_share'][B]:.3f}"),
        ("A|B spatial purity, observed", f"{ph['ab_purity']['observed']:.4f}"),
        ("A|B spatial purity, chance", f"{ph['ab_purity']['null_mean']:.4f}"),
        ("A|B spatial purity, z", f"{ph['ab_purity']['z_score']:.2f}"),
        (f"label agreement k_{K_BELOW:03d} vs unperturbed", f"{rec['label_agreement_below']:.4f}"),
        (f"ARI k_{K_BELOW:03d} vs unperturbed", f"{rec['ari_below']:.3f}"),
        (f"ARI k_{K_IDENTITY:03d} vs unperturbed", f"{rec['ari_identity']:.3f}"),
        ("SpatialAIM proposal, unperturbed run", ", ".join(f"{k}={v}" for k, v in sel["original"]["best"].items())),
        ("SpatialAIM proposal, perturbed run", ", ".join(f"{k}={v}" for k, v in sel["split"]["best"].items())),
    ]
    (out / "results.json").write_text(json.dumps(results, indent=1, default=float), encoding="utf-8")
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(out / "metrics.csv", index=False)
    logger.info("Wrote %s", out)
    print(pd.DataFrame(rows, columns=["metric", "value"]).to_string(index=False))
    print("\nscore delta (k_004 - k_003):",
          {k: round(v, 3) for k, v in results["score_delta_identity_vs_below"].items()})


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
