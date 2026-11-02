"""Parse GEO SOFT sample records into a tidy sample metadata table.

Group labels are taken from the ``treatment`` characteristic that the authors
submitted for each GSM record, not inferred from sample names. The sample-name
pattern is only used afterwards as an independent cross-check.
"""

import gzip
import re
from pathlib import Path

import pandas as pd

# GEO "treatment" characteristic -> (diet, microplastic) factor levels.
TREATMENT_TO_FACTORS = {
    "Chow Diet": ("CD", "Veh"),
    "Chow Diet + Microplastics": ("CD", "PE"),
    "MASH Diet": ("MD", "Veh"),
    "MASH Diet + Microplastics": ("MD", "PE"),
}


def read_soft_text(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return handle.read()


def split_sample_blocks(soft_text: str) -> dict[str, list[str]]:
    """Return {GSM accession: list of '!Sample_...' lines} for every sample entity."""
    blocks: dict[str, list[str]] = {}
    current = None
    for line in soft_text.splitlines():
        if line.startswith("^"):
            match = re.match(r"\^SAMPLE = (GSM\d+)", line)
            current = match.group(1) if match else None
            if current:
                blocks[current] = []
        elif current and line.startswith("!Sample_"):
            blocks[current].append(line)
    return blocks


def _field_values(lines: list[str], field: str) -> list[str]:
    prefix = f"!Sample_{field} = "
    return [line[len(prefix):].strip() for line in lines if line.startswith(prefix)]


def _key_value_fields(values: list[str]) -> dict[str, str]:
    """Turn ['tissue: Liver', 'batch: B1'] into {'tissue': 'Liver', 'batch': 'B1'}."""
    parsed = {}
    for value in values:
        key, _, rest = value.partition(":")
        parsed[key.strip()] = rest.strip()
    return parsed


def parse_sample_block(gsm: str, lines: list[str]) -> dict:
    """Extract the fields we use from one GSM record."""
    characteristics = _key_value_fields(_field_values(lines, "characteristics_ch1"))
    descriptions = _key_value_fields(_field_values(lines, "description"))
    relations = _key_value_fields(_field_values(lines, "relation"))

    counts_column = next(
        (v for k, v in descriptions.items() if k.startswith("column name in")), None
    )
    biosample_url = relations.get("BioSample", "")
    sra_url = relations.get("SRA", "")

    return {
        "geo_accession": gsm,
        "sample_title": _field_values(lines, "title")[0],
        "counts_column": counts_column,
        "library_name": descriptions.get("Library name"),
        "tissue": characteristics.get("tissue"),
        "genotype": characteristics.get("genotype"),
        "geo_treatment": characteristics.get("treatment"),
        "geo_batch": characteristics.get("batch"),
        "organism": _field_values(lines, "organism_ch1")[0],
        "molecule": _field_values(lines, "molecule_ch1")[0],
        "instrument": _field_values(lines, "instrument_model")[0],
        "library_strategy": _field_values(lines, "library_strategy")[0],
        "biosample": biosample_url.rsplit("/", 1)[-1] if biosample_url else None,
        "sra_experiment": sra_url.rsplit("=", 1)[-1] if sra_url else None,
        "submission_date": _field_values(lines, "submission_date")[0],
    }


def add_design_factors(metadata: pd.DataFrame) -> pd.DataFrame:
    """Add diet / microplastic / group columns from the GEO treatment field."""
    unknown = set(metadata["geo_treatment"]) - set(TREATMENT_TO_FACTORS)
    if unknown:
        raise ValueError(f"unrecognised GEO treatment values: {sorted(unknown)}")
    factors = metadata["geo_treatment"].map(TREATMENT_TO_FACTORS)
    out = metadata.copy()
    out["diet"] = [f[0] for f in factors]
    out["microplastic"] = [f[1] for f in factors]
    out["group"] = out["diet"] + "+" + out["microplastic"]
    out["replicate"] = out.groupby("group").cumcount() + 1
    return out


def group_from_column_name(column: str) -> str:
    """Infer the group from a counts column name such as 'MD_MP_3_S55'.

    Used only to cross-check the GEO treatment field, never as the primary label.
    """
    match = re.fullmatch(r"(CD|MD)(_MP)?_\d+_S\d+", column)
    if not match:
        raise ValueError(f"unexpected counts column name: {column}")
    return f"{match.group(1)}+{'PE' if match.group(2) else 'Veh'}"


def build_sample_metadata(soft_path: Path) -> pd.DataFrame:
    blocks = split_sample_blocks(read_soft_text(soft_path))
    records = [parse_sample_block(gsm, lines) for gsm, lines in blocks.items()]
    metadata = add_design_factors(pd.DataFrame.from_records(records))
    metadata["group_from_column_name"] = metadata["counts_column"].map(group_from_column_name)
    return metadata
