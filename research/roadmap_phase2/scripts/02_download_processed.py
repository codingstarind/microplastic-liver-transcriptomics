"""Step 02: download the deposited processed expression files for the
shortlisted series (no FASTQ). Rerunnable: existing files are skipped,
checksums are refreshed in metadata/download_manifest.csv.
"""
from p2_common import fetch

FILES = {
    "GSE316023": ["GSE316023_COUNT_NCOMMS-25-40823-T.csv.gz", "GSE316023_FPKM_NCOMMS-25-40823-T.csv.gz"],
    "GSE246328": ["GSE246328_GAN-DIO-NASH-Senescence_counts.txt.gz"],
    "GSE294629": ["GSE294629_processed_file.tsv.gz"],
    "GSE107787": ["GSE107787_RNAseq_expression_results.xlsx"],
    "GSE117134": ["GSE117134_casava_gene_expression.tsv.gz"],
    "GSE282998": ["GSE282998_GW7647_gene_counts.csv.gz"],
    # Track B status checks (small files only)
    "GSE245069": ["GSE245069_gene_expression.txt.gz"],
    "GSE277772": ["GSE277772_RAW.tar", "filelist.txt"],
    "GSE261893": ["GSE261893_readcount_Con_PE.txt.gz"],
    "GSE261725": ["GSE261725_readcount_Con_PS.txt.gz"],
    "GSE286116": ["GSE286116_readcount_Con_PP.txt.gz"],
    "GSE262399": ["GSE262399_readcount_Con_PB.txt.gz"],
}
MATRIX = {"GSE8295": "GSE8295_series_matrix.txt.gz"}  # normalized microarray values

if __name__ == "__main__":
    for acc, names in FILES.items():
        for n in names:
            fetch(acc, "suppl/", n)
    for acc, n in MATRIX.items():
        fetch(acc, "matrix/", n)
    fetch("GSE225695", "suppl/", "GSE225695_mRNA_Expression_Profiling.xlsx")
    # MGI marker list (symbols + synonyms) for renamed-gene rescue; ~86 MB, dated snapshot
    from p2_common import download, ds_dir, record_manifest
    mgi_url = "https://www.informatics.jax.org/downloads/reports/MRK_List1.rpt"
    dest = ds_dir("MGI") / "MRK_List1.rpt"
    download(mgi_url, dest)
    record_manifest("MGI", mgi_url, dest)
