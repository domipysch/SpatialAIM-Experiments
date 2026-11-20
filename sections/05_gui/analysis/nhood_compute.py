import numpy as np, pandas as pd, anndata as ad, squidpy as sq
from spatialaim_experiments import paths

D = paths.GUI
OUT = paths.GUI_ANALYSIS / "nhood_k016"
st = ad.read_h5ad(D / "ST/11_human-crc-p2.h5ad", backed="r")
m = ad.read_h5ad(D / "Result/leiden/nearest_centroid/k_016/spot_to_state_mapping_soft.h5ad")
assert (st.obs_names == m.obs_names).all()
X = m.X if isinstance(m.X, np.ndarray) else m.X.toarray()
a = ad.AnnData(obs=pd.DataFrame(index=m.obs_names))
a.obsm["spatial"] = np.asarray(st.obsm["spatial"])
hard = X.argmax(1)
a.obs["state"] = pd.Categorical(hard.astype(str), categories=[str(i) for i in range(X.shape[1])])
print("counts", np.bincount(hard, minlength=16).tolist())
sq.gr.spatial_neighbors_knn(a, n_neighs=6)
sq.gr.nhood_enrichment(a, cluster_key="state", n_perms=500, seed=0, n_jobs=1, show_progress_bar=False)
r = a.uns["state_nhood_enrichment"]
cats = list(a.obs["state"].cat.categories)
z = pd.DataFrame(r["zscore"], index=cats, columns=cats)
c = pd.DataFrame(r["count"], index=cats, columns=cats)
OUT.mkdir(parents=True, exist_ok=True)
z.to_csv(OUT / "nhood_z.csv"); c.to_csv(OUT / "nhood_count.csv")
print(z.round(0).astype(int).to_string())
