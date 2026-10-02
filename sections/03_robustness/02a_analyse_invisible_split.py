import argparse
import json
import logging
import anndata as ad
import numpy as np
import pandas as pd
from spatialaim_experiments import paths
from spatialaim_experiments.style import MAPPERS

logger = logging.getLogger(__name__)


RUN_SPLIT = "invisible_split"  # run folder under Data/02_Results/robustness

A, B = "Glutamatergic_A", "Glutamatergic_B"
K_IDENTITY = 4  # the split annotation's own four types


def state_of(run_dir, k: int, type_name: str) -> int:
    """The state the annotated type ``type_name`` belongs to at the cut ``k``."""
    obs = ad.read_h5ad(run_dir / "start_clustering.h5ad").obs
    start = int(obs.loc[obs["start_cluster_name"].astype(str) == type_name, "start_cluster"].iloc[0])
    cut = pd.read_csv(run_dir / f"k_{k:03d}" / "start_cluster_to_state.csv")
    return int(cut.loc[cut["start_cluster"] == start, "state"].iloc[0])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mapper", required=True, choices=MAPPERS,
                    help="which mapper's runs to read")
    mapper = ap.parse_args().mapper
    run_dir = paths.ROBUSTNESS / RUN_SPLIT / mapper

    # the phantom division at the identity cut, as a hard spot assignment
    P = np.asarray(ad.read_h5ad(run_dir / f"k_{K_IDENTITY:03d}" /
                                "spot_to_state_mapping_soft.h5ad").X, dtype=float)
    hard = np.argmax(P, axis=1)
    sa, sb = state_of(run_dir, K_IDENTITY, A), state_of(run_dir, K_IDENTITY, B)
    spot_share = {A: float((hard == sa).mean()), B: float((hard == sb).mean())}

    selection = json.loads((run_dir / "k_selection.json").read_text(encoding="utf-8"))
    results = {
        "mapper": mapper,
        "spot_share": spot_share,
        "k_proposed": {name: (None if hit is None else int(hit["k"]))
                       for name, hit in selection["best"].items()},
        "scores": pd.read_csv(run_dir / "k_comparison.csv").to_dict(orient="list"),
    }

    out = paths.ROBUSTNESS / "_analysis" / mapper
    out.mkdir(parents=True, exist_ok=True)
    (out / "results.json").write_text(json.dumps(results, indent=1, default=float), encoding="utf-8")
    logger.info("Wrote %s", out / "results.json")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
