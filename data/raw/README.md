# data/raw

Files here are exactly as downloaded and must never be edited. They are not
committed to git; recreate them with:

```bash
python scripts/01_download_data.py
```

The script verifies each file against the SHA-256 recorded in
`metadata/data_provenance.json` and refuses to overwrite a file that differs.

| File | Source | What it is |
|---|---|---|
| `GSE294180_MP_raw_counts_gene_level.csv.gz` | GEO GSE294180 supplementary | Gene-level counts from Salmon + tximport (19 samples) |
| `GSE294180_family.soft.gz` | GEO GSE294180 SOFT | Series, platform and sample records (metadata source) |
| `gencode.vM35.annotation.gtf.gz` | GENCODE mouse M35 (EBI) | Gene annotation used to add symbols / biotypes |
