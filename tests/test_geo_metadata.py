import pandas as pd
import pytest

from mpliver.geo_metadata import (
    add_design_factors,
    group_from_column_name,
    parse_sample_block,
    split_sample_blocks,
)

SOFT_SNIPPET = """^SERIES = GSE294180
!Series_title = test
^SAMPLE = GSM1
!Sample_title = Liver_CD_MP_1
!Sample_submission_date = Apr 09 2025
!Sample_organism_ch1 = Mus musculus
!Sample_characteristics_ch1 = tissue: Liver
!Sample_characteristics_ch1 = genotype: WT
!Sample_characteristics_ch1 = treatment: Chow Diet + Microplastics
!Sample_characteristics_ch1 = batch: B1
!Sample_molecule_ch1 = total RNA
!Sample_description = Library name: S48
!Sample_description = column name in MP_raw_counts_gene_level.csv: CD_MP_1_S48
!Sample_instrument_model = NextSeq 2000
!Sample_library_strategy = RNA-Seq
!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN47862335
!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX28316101
^SAMPLE = GSM2
!Sample_title = Liver_MD_1
"""


def test_split_sample_blocks_finds_each_sample():
    blocks = split_sample_blocks(SOFT_SNIPPET)
    assert list(blocks) == ["GSM1", "GSM2"]
    assert not any(line.startswith("!Series") for line in blocks["GSM1"])


def test_parse_sample_block_extracts_fields():
    record = parse_sample_block("GSM1", split_sample_blocks(SOFT_SNIPPET)["GSM1"])
    assert record["counts_column"] == "CD_MP_1_S48"
    assert record["geo_treatment"] == "Chow Diet + Microplastics"
    assert record["geo_batch"] == "B1"
    assert record["biosample"] == "SAMN47862335"
    assert record["sra_experiment"] == "SRX28316101"


def test_add_design_factors_maps_treatments():
    meta = pd.DataFrame({"geo_treatment": ["Chow Diet", "MASH Diet + Microplastics",
                                           "MASH Diet + Microplastics"]})
    out = add_design_factors(meta)
    assert out["group"].tolist() == ["CD+Veh", "MD+PE", "MD+PE"]
    assert out["replicate"].tolist() == [1, 1, 2]


def test_add_design_factors_rejects_unknown_treatment():
    with pytest.raises(ValueError):
        add_design_factors(pd.DataFrame({"geo_treatment": ["High fat"]}))


@pytest.mark.parametrize("column, group", [
    ("CD_1_S39", "CD+Veh"),
    ("CD_MP_5_S52", "CD+PE"),
    ("MD_2_S44", "MD+Veh"),
    ("MD_MP_3_S55", "MD+PE"),
])
def test_group_from_column_name(column, group):
    assert group_from_column_name(column) == group


def test_group_from_column_name_rejects_unexpected():
    with pytest.raises(ValueError):
        group_from_column_name("Sample_7")
