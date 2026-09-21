# Supplement

```bash
python supplement/build_dataset_tables.py
```

Builds `SpatialAIM-Overleaf/supplement/supplementary_datasets.xlsx` — three sheets
(scRNA references, spatial slices, pairs) from the two dataset index CSVs,
restricted to the 92 whole pairs the paper uses.

The `Citation` column is left empty on purpose: the indices carry a URL and a
first author but no bib key, so the keys are filled in by hand.
