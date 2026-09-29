#!/usr/bin/env python3
"""Export eleven primary CSV views without changing saved numerical tokens."""
import argparse
import csv
import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path


REGIONS = {"LC": "LC", "SNVTA": "SN–VTA"}
LABELS = {
    **REGIONS,
    "lme_intercept": "Intercept",
    "lme_slope": "Slope",
    "mci_to_dementia": "MCI-to-dementia",
    "any_worsening": "Any progression",
    "primary_individual_index": "Primary separate LC and SN–VTA models",
    "same_frame_base": "Without TPV adjustment",
    "same_frame_tpv": "With TPV adjustment",
    "binary": "Binary amyloid",
    "centiloid": "Continuous amyloid (Centiloid)",
    "(Intercept)": "Intercept",
    "age_baseline_c": "Baseline age (years)",
    "mri_time_years_c": "Imaging time (years)",
    "sexM": "Male sex (reference: female)",
    "diagnosis_baselineMCI": "Baseline MCI (reference: CN)",
    "diagnosis_baselineAD": "Baseline AD (reference: CN)",
    "hemisphere_c": "Right versus left hemisphere",
    "field_strength3T": "3 T (reference: 1.5 T)",
    "education_years_c": "Education (years)",
    "apoe4_count": "APOE ε4 allele count",
    "efc_c": "Entropy focus criterion",
}
LABEL_FIELDS = {
    "nucleus", "focal_nucleus", "branch", "model", "term", "endpoint",
    "outcome", "context", "model_variant", "encoding", "effect_role",
}
PRIMARY_MODELS = {
    "concurrent": (
        "original", "not_applicable", "Concurrent (hemisphere-specific CR outcome)",
    ),
    "prospective": ("primary", "180_1095", "Prospective (180–1,095 days)"),
    "lag": ("full_lag", "not_applicable", "All-assessment lag model"),
}


def region_label(value):
    try:
        return REGIONS[value]
    except KeyError as error:
        raise ValueError("Unknown CR Region: " + value) from error


def label(value):
    if value.startswith("echo_fTE_"):
        return "Echo time " + value[len("echo_fTE_"):] + " s versus recorded reference"
    return LABELS.get(value, value)


def public_label(field, value, row):
    if field in ("nucleus", "focal_nucleus"):
        return region_label(value)
    if field == "effect_role":
        cr_region = region_label(row["nucleus"])
        return {
            "local_pathology_time": cr_region + " CR × amyloid × time",
            "plsc_pathology_time": cr_region + "-related brain pattern × amyloid × time",
            "tpv_pathology_time": "TPV × amyloid × time",
        }[value]
    return label(value)


def columns(*items):
    return [(item, item) if isinstance(item, str) else item for item in items]


def primary_f3(row, stage):
    variant, window, _ = PRIMARY_MODELS[stage]
    return (
        row["model_type"] == stage
        and row["variant"] == variant
        and row["window"] == window
        and row["comparison_model"] == "single"
        and row["pairing"] == "primary"
    )


def model_type_label(row, stage):
    if not primary_f3(row, stage):
        raise ValueError("Expected the native primary Figure 3 identity: " + stage)
    return PRIMARY_MODELS[stage][2]


# filename -> native path, stage, exact source/destination columns, row filter.
VIEWS = {}


def add(name, source, stage, cols, select=lambda row: True):
    VIEWS[name + ".csv"] = (source, stage, cols, select)


add(
    "icc",
    "figure1/icc.csv",
    "icc",
    columns(
        ("nucleus", "Region"),
        ("field_strength", "Field_strength"),
        ("hemisphere", "Hemisphere"),
        ("n_participants", "N"),
        ("n_observations", "Observations"),
        ("icc", "ICC"),
        ("ci_low", "CI_lower"),
        ("ci_high", "CI_upper"),
    ),
    lambda r: r["variant"] == "primary",
)
add(
    "hypopigmentation_associations",
    "figure1/gross_coefficients.csv",
    "pathology",
    columns(
        ("nucleus", "Region"),
        ("rating", "Rating"),
        ("n", "N"),
        ("n_absent", "Rating_absent_N"),
        ("n_present", "Rating_present_N"),
        "beta",
        ("se", "SE"),
        ("ci_low", "CI_lower"),
        ("ci_high", "CI_upper"),
        "P",
    ),
    lambda r: r["focal"].upper() == "TRUE",
)
add(
    "demographic_associations",
    "figure2/coefficients.csv",
    "figure2",
    columns(
        ("nucleus", "Region"),
        ("term", "Variables"),
        "beta",
        ("se", "SE"),
        ("statistic", "t"),
        "P",
    ),
    lambda r: r["variant"] == "primary" and r["record_type"] == "fixed_coefficient",
)
add(
    "plsc_results",
    "figure4/plsc_global_results.csv",
    "plsc",
    columns(
        ("branch", "Model type"),
        ("model", "Parameter"),
        ("n_subjects", "N"),
        ("brain_behavior_r", "r"),
        "P",
    ),
    lambda r: r["branch"] in ("LC", "SNVTA")
    and r["model"] in ("lme_intercept", "lme_slope"),
)
add(
    "molecular_pet_correspondence",
    "figure4/molecular_results.csv",
    "molecular",
    columns(
        ("branch", "Model type"),
        ("target", "PET_target"),
        ("tracer", "PET_tracer"),
        ("source", "PET_source"),
        ("pet_analysis_family", "Family"),
        ("rho", "Spearman_rho"),
        "P",
        "P_fdr",
    ),
    lambda r: r["branch"] in ("LC", "SNVTA") and r["model"] == "lme_intercept",
)
add(
    "clinical_progression",
    "figure5/progression_results.csv",
    "progression",
    columns(
        ("endpoint", "Outcome"),
        ("context", "Model type"),
        ("n_subjects", "N"),
        ("n_events", "Events"),
        ("focal_nucleus", "Region"),
        "beta",
        ("se", "SE"),
        ("odds_ratio", "OR"),
        ("ci_low", "CI_lower"),
        ("ci_high", "CI_upper"),
        "P",
    ),
    lambda r: r["context"] == "primary_individual_index"
    and r["nucleus"] in ("LC", "SNVTA"),
)
add(
    "mediation_effects",
    "figure5/mediation_effects.csv",
    "mediation",
    columns(
        ("nucleus", "Region"),
        ("outcome", "Outcome"),
        ("model_variant", "TPV_adjustment"),
        ("n", "N"),
        ("effect", "Effect"),
        ("estimate", "beta"),
        ("se", "SE"),
        ("ci_low", "CI_lower"),
        ("ci_high", "CI_upper"),
        "P",
    ),
    lambda r: r["definition"] == "primary_trajectory"
    and r["model_variant"] in ("same_frame_base", "same_frame_tpv")
    and r["effect"] in ("indirect", "direct", "total"),
)
add(
    "amyloid_memory_interactions",
    "figure5/memory_focal_effects.csv",
    "amyloid",
    columns(
        ("nucleus", "Region"),
        ("encoding", "Amyloid_encoding"),
        ("n_subjects", "N"),
        ("n_rows", "Observations"),
        ("effect_role", "Interaction"),
        ("estimate", "beta"),
        ("se", "SE"),
        ("ci_low", "CI_lower"),
        ("ci_high", "CI_upper"),
        "P",
    ),
    lambda r: r["scope"] == "primary"
    and r["model_variant"] == "R2_local_plsc_tpv"
    and r["effect_role"]
    in ("local_pathology_time", "plsc_pathology_time", "tpv_pathology_time"),
)

FOCAL_COLUMNS = columns(
    ("variable_name", "Outcome"),
    ("nucleus", "Region"),
    ("variant", "Model type"),
    ("term", "Term"),
    "beta",
    "SE",
    ("ci_low", "CI_lower"),
    ("ci_high", "CI_upper"),
    ("statistic", "t"),
    "df",
    "P",
    "P_fdr",
    ("n_participants", "N"),
    ("n_rows", "Observations"),
)
LAG_COLUMNS = columns(
    ("variable_name", "Outcome"),
    ("nucleus", "Region"),
    ("variant", "Model type"),
    ("n_participants", "N"),
    ("n_scans", "MRI scans"),
    ("n_rows", "Observations"),
    ("n_intervals", "Interval number"),
    ("significant_intervals", "Direction"),
    ("significant_intervals", "Start (years)"),
    ("significant_intervals", "End (years)"),
)
FOCAL_HEADERS = [dest for _, dest in FOCAL_COLUMNS]
LAG_HEADERS = [dest for _, dest in LAG_COLUMNS]
NUMBER = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
INTERVAL = re.compile(
    r"(positive|negative): (" + NUMBER + r") to (" + NUMBER + r") years"
)
MISSING = ("", "NA", "NaN", "nan")


def parse_intervals(count, text):
    """Copy interval tokens; count is a per-curve total, never a row ordinal."""
    if count in MISSING:
        if text not in MISSING and text not in ("Not estimated", "Unavailable"):
            raise ValueError(
                "Unestimated interval count has an estimated interval value"
            )
        return [("NA", "NA", "NA")]
    if not re.fullmatch(r"[0-9]+", count):
        raise ValueError("Invalid interval count: " + repr(count))
    number = int(count)
    if number == 0:
        if text != "None":
            raise ValueError("Zero intervals must have the saved None marker")
        return [("None", "NA", "NA")]
    parts = text.split("; ")
    if len(parts) != number:
        raise ValueError("Declared interval count does not match saved intervals")
    values = []
    for part in parts:
        match = INTERVAL.fullmatch(part)
        if not match:
            raise ValueError("Malformed saved interval: " + repr(part))
        direction, start, end = match.groups()
        if Decimal(start) > Decimal(end):
            raise ValueError("Interval start exceeds end")
        values.append((direction, start, end))
    if len(set(values)) != len(values):
        raise ValueError("Duplicated saved interval")
    return values


def clean_token(value):
    return "NA" if value in MISSING else value


def figure3_rows(rows, stage):
    """Copy native primary associations or expand their saved lag intervals."""
    headers = LAG_HEADERS if stage == "lag" else FOCAL_HEADERS
    output, keys = [], set()
    for row in rows:
        cr_region = region_label(row["nucleus"])
        identity = {
            "Outcome": row["variable_name"],
            "Region": cr_region,
            "Model type": model_type_label(row, stage),
        }
        key = tuple(identity.values())
        if key in keys:
            raise ValueError("Duplicate Figure 3 scientific key: " + repr(key))
        keys.add(key)
        if stage == "lag":
            result = {dest: clean_token(row[src]) for src, dest in LAG_COLUMNS[:7]}
            result.update(identity)
            intervals = parse_intervals(row["n_intervals"], row["significant_intervals"])
            for direction, start, end in intervals:
                output.append(dict(result, **{
                    "Direction": direction, "Start (years)": start, "End (years)": end,
                }))
        else:
            result = {dest: clean_token(row[src]) for src, dest in FOCAL_COLUMNS}
            result.update(identity, Term=cr_region + " CR")
            output.append(result)
    return headers, [{field: row[field] for field in headers} for row in output]


for stage in ("concurrent", "prospective"):
    add(
        "Figure3_" + stage + "_associations",
        "figure3/" + stage + "_associations.csv",
        stage,
        FOCAL_COLUMNS,
        lambda r, stage=stage: primary_f3(r, stage),
    )
add(
    "Figure3_lag_interval_summary",
    "figure3/lag_interval_summary.csv",
    "lag",
    LAG_COLUMNS,
    lambda r: primary_f3(r, "lag"),
)

KEYS = {
    "icc.csv": ("nucleus", "field_strength", "hemisphere"),
    "hypopigmentation_associations.csv": ("nucleus", "rating"),
    "demographic_associations.csv": ("nucleus", "term"),
    "plsc_results.csv": ("branch", "model"),
    "molecular_pet_correspondence.csv": ("branch", "pet_id"),
    "clinical_progression.csv": ("endpoint", "nucleus"),
    "mediation_effects.csv": ("nucleus", "outcome", "model_variant", "effect"),
    "amyloid_memory_interactions.csv": ("nucleus", "encoding", "effect_role"),
}
for name in VIEWS:
    if name.startswith("Figure3_"):
        KEYS[name] = ("outcome", "nucleus", "variant", "window", "comparison_model")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(
            reader.fieldnames
        ):
            raise ValueError("Missing or duplicate CSV header: " + str(path))
        rows = list(reader)
        if any(
            None in row or any(value is None for value in row.values()) for row in rows
        ):
            raise ValueError("Ragged CSV: " + str(path))
        return reader.fieldnames, rows


def export_views(output, source_root, stages=None, receipt=None):
    """Validate native sources, then write immutable aggregate views atomically."""
    output = Path(output).resolve()
    source_root = Path(source_root).resolve(strict=True)
    selected = {
        name: value
        for name, value in VIEWS.items()
        if stages is None or value[1] in stages
    }
    if not selected:
        return {
            "schema": "main_csv_views_v1",
            "files": {},
            "scientific_calculations": 0,
        }
    if output.exists() and (
        any(p.is_dir() for p in output.iterdir())
        or any(p.name not in VIEWS for p in output.iterdir())
    ):
        raise ValueError("Ordinary result directory contains unexpected files")
    sources, results = {}, {}
    # Validate every source and construct every view before writing any output.
    for name, (source, stage, cols, select) in selected.items():
        path = (source_root / source).resolve(strict=True)
        if output == path.parent or output in path.parents:
            raise ValueError("Output cannot contain a selected source")
        digest = sha(path)
        fields, rows = read_csv(path)
        absent = {src for src, _ in cols} - set(fields)
        if absent:
            raise ValueError(
                "Missing source columns in " + source + ": " + ", ".join(sorted(absent))
            )
        if name.startswith("Figure3_"):
            for row in rows:
                model_type_label(row, stage)
        kept = [(index, row) for index, row in enumerate(rows, 1) if select(row)]
        keys = [tuple(row[field] for field in KEYS[name]) for _, row in kept]
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate scientific key in selected view: " + name)
        if name.startswith("Figure3_"):
            fields, transformed = figure3_rows([row for _, row in kept], stage)
            output_source_rows = [
                index
                for index, row in kept
                for _ in (
                    parse_intervals(row["n_intervals"], row["significant_intervals"])
                    if stage == "lag" else [None]
                )
            ]
        else:
            fields = [dest for _, dest in cols]
            output_source_rows = [index for index, _ in kept]
            transformed = []
            for _, row in kept:
                result = {}
                for src, dest in cols:
                    value = clean_token(row[src])
                    if name == "hypopigmentation_associations.csv" and src == "rating":
                        value = {
                            "LC": "LC hypopigmentation",
                            "SNVTA": "SN hypopigmentation",
                            "SN": "SN hypopigmentation",
                        }[value]
                    elif src in LABEL_FIELDS and value != "NA":
                        value = public_label(src, value, row)
                    result[dest] = value
                transformed.append(result)
        results[name] = (fields, transformed)
        sources[name] = {
            "native_file": source,
            "path": str(path),
            "sha256": digest,
            "source_rows": len(rows),
            "selected_source_rows": [i for i, _ in kept],
            "output_source_rows": output_source_rows,
            "rows": len(transformed),
            "columns": fields,
            "column_sources": dict((dest, src) for src, dest in cols),
        }
    output.mkdir(parents=True, exist_ok=True)
    for name, (fields, rows) in results.items():
        target = output / name
        temporary = output / (name + ".tmp")
        with temporary.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        if target.exists():
            unchanged = sha(target) == sha(temporary)
            temporary.unlink()
            if not unchanged:
                raise ValueError(
                    "Existing view differs; use a fresh output root: " + name
                )
        else:
            temporary.replace(target)
        sources[name]["output_sha256"] = sha(target)
    record = {
        "schema": "main_csv_views_v1",
        "exporter_sha256": sha(Path(__file__)),
        "scientific_calculations": 0,
        "files": sources,
    }
    if receipt:
        receipt = Path(receipt).resolve()
        if receipt.parent == output or output in receipt.parents:
            raise ValueError("Technical receipt must be outside ordinary results")
        receipt.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(record, indent=2, sort_keys=True) + "\n"
        if not receipt.exists() or receipt.read_text() != text:
            receipt.write_text(text)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    record = export_views(
        args.out, args.source_root, receipt=args.receipt
    )
    print(
        "Exported",
        len(record["files"]),
        "aggregate CSV views; zero scientific calculations",
    )


if __name__ == "__main__":
    main()
