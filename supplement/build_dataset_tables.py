"""Supplementary dataset tables -> one xlsx with three sheets.

Restricted to the 92 whole pairs the paper actually uses
(``pairs_only_raw_counts_11_12_whole.csv``), so the 244-pair database and the
tiled 11_<n>/12_<n> variants stay out.

Sheets
------
scRNA   one row per reference used (10)
ST      one row per spatial slice used (92)
Pairs   one row per pair, the two names plus the shared-gene count

Everything comes from the two index CSVs; the ``Citation`` column is left empty
for the bib keys to be filled in by hand, since the indices carry a URL and a
first author but no reference key.

    python supplement/build_dataset_tables.py
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from aim_experiments import paths

logger = logging.getLogger(__name__)

DATASETS = paths.DATASETS
PAIRS_CSV = paths.PAIRS_CSV
OUT = paths.PAPER / "supplement" / "supplementary_datasets.xlsx"

# index column -> column heading in the supplement. Renamed rather than kept as
# they are: "Quelle" and "Autor" are German, and "Links" is the download
# location, which is a different thing from the publication.
SC_COLUMNS = {
    "Name": "Dataset",
    "Tissue": "Tissue",
    "Autor": "Source publication",
    "Quelle": "Publication URL",
    "Citation": "Citation",
    "Tech": "Technology",
    "CellCount": "Cells",
    "GeneCount": "Genes",
    "CellTypeKey0": "Annotation 1",
    "NumberCellTypes0": "Cell types 1",
    "CellTypeKey1": "Annotation 2",
    "NumberCellTypes1": "Cell types 2",
    "Links": "Data availability",
}
ST_COLUMNS = {
    "Name": "Slice",
    "Reference": "Paired reference",
    "Tissue": "Tissue",
    "Autor": "Source publication",
    "Quelle": "Publication URL",
    "Citation": "Citation",
    "Tech": "Technology",
    "SpotCount": "Cells",
    "GeneCount": "Genes",
    "Links": "Data availability",
}
PAIR_COLUMNS = {
    "PairID": "Pair",
    "scName": "scRNA-seq reference",
    "stName": "Spatial slice",
    "NumberSharedGenes": "Shared genes",
}


def _select(index: pd.DataFrame, names, columns: dict) -> pd.DataFrame:
    table = index[index["Name"].isin(names)].copy()
    table["Citation"] = ""
    missing = [c for c in columns if c not in table.columns]
    assert not missing, f"index is missing {missing}"
    return table[list(columns)].rename(columns=columns).sort_values(
        by=list(columns.values())[0]
    )


def build() -> dict[str, pd.DataFrame]:
    pairs = pd.read_csv(PAIRS_CSV)
    sc_index = pd.read_csv(DATASETS / "scRNA" / "index.csv")
    st_index = pd.read_csv(DATASETS / "ST" / "index.csv")

    sc = _select(sc_index, set(pairs["scName"]), SC_COLUMNS)

    # every slice is paired with exactly one reference in this corpus, so the
    # mapping is well defined and worth carrying: the ST index does not record it
    st_index = st_index.merge(
        pairs[["stName", "scName"]].rename(
            columns={"stName": "Name", "scName": "Reference"}
        ),
        on="Name",
        how="left",
    )
    st = _select(st_index, set(pairs["stName"]), ST_COLUMNS)

    pair_table = pairs[list(PAIR_COLUMNS)].rename(columns=PAIR_COLUMNS)
    return {"scRNA": sc, "ST": st, "Pairs": pair_table}


def write(sheets: dict[str, pd.DataFrame], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, table in sheets.items():
            table.to_excel(writer, sheet_name=name, index=False)
            worksheet = writer.sheets[name]
            worksheet.freeze_panes = "A2"
            for i, column in enumerate(table.columns, start=1):
                width = max(len(str(column)), *(len(str(v)) for v in table[column]))
                worksheet.column_dimensions[
                    worksheet.cell(row=1, column=i).column_letter
                ].width = min(width + 2, 60)


def main() -> None:
    sheets = build()
    write(sheets, OUT)
    for name, table in sheets.items():
        logger.info("%-6s %3d rows x %d columns", name, len(table), table.shape[1])
    logger.info("-> %s", OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
