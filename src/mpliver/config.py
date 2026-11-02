"""Project-wide paths and dataset constants.

Every script imports paths from here so that nothing depends on the current
working directory.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = REPO_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
METADATA_DIR = REPO_ROOT / "metadata"
FIGURES_DIR = REPO_ROOT / "figures"
TABLES_DIR = REPO_ROOT / "results" / "tables"
REPORTS_DIR = REPO_ROOT / "reports"

GEO_ACCESSION = "GSE294180"
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294180"

# Files deposited by the authors (the only processed data GEO provides) plus the
# machine-readable series record that holds the per-sample metadata.
SOURCE_FILES = {
    "counts": {
        "url": f"{GEO_FTP_BASE}/suppl/GSE294180_MP_raw_counts_gene_level.csv.gz",
        "filename": "GSE294180_MP_raw_counts_gene_level.csv.gz",
        "description": "Gene-level counts (Salmon -> tximport) deposited by the authors",
    },
    "family_soft": {
        "url": f"{GEO_FTP_BASE}/soft/GSE294180_family.soft.gz",
        "filename": "GSE294180_family.soft.gz",
        "description": "GEO SOFT family file: series, platform and all sample records",
    },
    # Same reference release the authors quantified against (GEO: "GRCm39 M35").
    # Used only to attach gene symbols, biotypes and chromosomes to Ensembl IDs.
    "gencode_gtf": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/"
               "gencode.vM35.annotation.gtf.gz",
        "filename": "gencode.vM35.annotation.gtf.gz",
        "description": "GENCODE mouse release M35 comprehensive gene annotation (GRCm39)",
    },
}

COUNTS_FILE = RAW_DIR / SOURCE_FILES["counts"]["filename"]
SOFT_FILE = RAW_DIR / SOURCE_FILES["family_soft"]["filename"]
GTF_FILE = RAW_DIR / SOURCE_FILES["gencode_gtf"]["filename"]
GENE_ANNOTATION_FILE = METADATA_DIR / "gene_annotation_gencode_vM35.csv.gz"
PROVENANCE_FILE = METADATA_DIR / "data_provenance.json"
SAMPLE_METADATA_FILE = METADATA_DIR / "sample_metadata.csv"

# Order and display labels for the 2 x 2 factorial design.
# Colours: categorical slots 1-4 of a palette validated for colour-vision
# deficiency. Two of them are light on white, so markers always get a dark edge,
# and microplastic status is also encoded by marker shape (never colour alone).
GROUP_ORDER = ["CD+Veh", "CD+PE", "MD+Veh", "MD+PE"]
GROUP_COLORS = {
    "CD+Veh": "#2a78d6",
    "CD+PE": "#eb6834",
    "MD+Veh": "#1baf7a",
    "MD+PE": "#eda100",
}
GROUP_MARKERS = {"CD+Veh": "o", "CD+PE": "^", "MD+Veh": "o", "MD+PE": "^"}
GROUP_LABELS = {
    "CD+Veh": "Control diet + vehicle",
    "CD+PE": "Control diet + PE",
    "MD+Veh": "MASH diet + vehicle",
    "MD+PE": "MASH diet + PE",
}
