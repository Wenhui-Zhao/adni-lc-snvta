#!/usr/bin/env python3
"""Compare direct-intercept PLSC saliences with cortical PET maps.

Use symmetric ten-neighbour Gaussian weights, singleton Moran randomization,
average-rank Spearman correlations and pooled-draw empirical probabilities.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
from typing import Any, Optional, Union
import numpy as np
import pandas as pd
import scipy
from scipy.linalg import eigh, helmert
from scipy.sparse.csgraph import connected_components
from scipy.stats import rankdata


def build_spatial_weights(
    distance: np.ndarray, config: dict[str, Any]
) -> tuple[np.ndarray, dict[str, Any]]:
    distance = np.asarray(distance, dtype=float)
    n = distance.shape[0]
    if distance.shape != (n, n) or n < 4:
        raise ValueError(f"Distance matrix has invalid shape: {distance.shape}")
    k = int(config["atlas"].get("weight_k_neighbors", min(10, n - 1)))
    k = max(1, min(k, n - 1))

    order = np.argsort(distance, axis=1)
    neighbor_mask = np.zeros((n, n), dtype=bool)
    for i in range(n):
        nbr = [j for j in order[i] if j != i][:k]
        neighbor_mask[i, nbr] = True
    neighbor_mask |= neighbor_mask.T
    positive_neighbor_distances = distance[neighbor_mask & (distance > 0)]
    if positive_neighbor_distances.size == 0:
        raise ValueError("No positive neighbor distances are available")

    bandwidth = float(np.median(positive_neighbor_distances))
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError(f"Invalid spatial-weight bandwidth: {bandwidth}")

    weights = np.zeros_like(distance, dtype=float)
    weights[neighbor_mask] = np.exp(-0.5 * (distance[neighbor_mask] / bandwidth) ** 2)
    weights = 0.5 * (weights + weights.T)
    np.fill_diagonal(weights, 0.0)
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError("Spatial-weight matrix is invalid")
    n_components, labels = connected_components((weights > 0).astype(int), directed=False)
    if n_components > int(config["atlas"].get("max_weight_components", 1)):
        counts = np.bincount(labels).tolist()
        raise ValueError(f"Spatial-weight graph has {n_components} components with sizes {counts}")
    return weights, {
        "weight_kernel": "gaussian",
        "weight_k_neighbors": k,
        "weight_bandwidth": bandwidth,
        "weight_n_components": int(n_components),
        "weight_sum": float(weights.sum()),
        "weight_hash": stable_hash(np.round(weights, 10).tolist(), length=20),
    }


def moran_i(values: np.ndarray, weights: np.ndarray) -> float:
    x = np.asarray(values, dtype=float).reshape(-1)
    w = np.asarray(weights, dtype=float)
    if w.shape != (x.size, x.size):
        raise ValueError("Value and spatial-weight dimensions differ")
    finite = np.isfinite(x)
    if not finite.all():
        raise ValueError("Moran's I requires finite values")
    z = x - x.mean()
    denom = float(z @ z)
    s0 = float(w.sum())
    if denom <= 0 or s0 <= 0:
        return float("nan")
    return float((x.size / s0) * ((z @ w @ z) / denom))


@dataclass(frozen=True)
class MoranBasis:
    eigenvectors: np.ndarray
    eigenvalues: np.ndarray
    weights: np.ndarray
    reconstruction_rmse: float


def compute_moran_basis(weights: np.ndarray) -> MoranBasis:
    """Return the complete orthonormal Moran basis for the mean-zero subspace."""
    w = np.asarray(weights, dtype=np.float64)
    if w.ndim != 2 or w.shape[0] != w.shape[1]:
        raise ValueError("Spatial weights must be square")
    n = w.shape[0]
    if n < 4:
        raise ValueError("At least four spatial observations are required")
    if not np.isfinite(w).all():
        raise ValueError("Spatial weights contain non-finite values")
    if not np.allclose(w, w.T, atol=1e-10, rtol=1e-8):
        raise ValueError("Spatial weights must be symmetric")
    w = 0.5 * (w + w.T)
    np.fill_diagonal(w, 0.0)

    # The transposed Helmert matrix spans the centered space.
    q = helmert(n, full=False).T
    restricted = q.T @ w @ q
    evals, evecs = eigh(restricted, check_finite=True)
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    mem = q @ evecs[:, order]
    orth_error = np.max(np.abs(mem.T @ mem - np.eye(n - 1)))
    center_error = np.max(np.abs(mem.T @ np.ones(n)))
    if orth_error > 1e-7 or center_error > 1e-7:
        raise RuntimeError(
            f"Failed to construct a stable Moran basis: orth_error={orth_error}, center_error={center_error}"
        )
    # Verify reconstruction of a centered vector.
    probe = np.arange(n, dtype=float) - (n - 1) / 2
    recon = mem @ (mem.T @ probe)
    rmse = float(np.sqrt(np.mean((recon - probe) ** 2)))
    return MoranBasis(eigenvectors=mem, eigenvalues=evals, weights=w, reconstruction_rmse=rmse)


def moran_randomize(
    values: np.ndarray,
    basis: MoranBasis,
    *,
    n_rep: int,
    random_state: Optional[Union[int, np.random.Generator]] = None,
) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64).reshape(-1)
    mem = basis.eigenvectors
    if x.size != mem.shape[0]:
        raise ValueError("Map length and Moran basis size differ")
    if not np.isfinite(x).all():
        raise ValueError("Moran randomization requires finite map values")
    if n_rep <= 0:
        raise ValueError("n_rep must be positive")
    if np.nanstd(x, ddof=1) <= 0:
        raise ValueError("Cannot randomize a constant map")
    rng = (
        random_state
        if isinstance(random_state, np.random.Generator)
        else np.random.default_rng(random_state)
    )
    centered = x - x.mean()
    coefficients = mem.T @ centered
    signs = rng.choice(np.array([-1.0, 1.0]), size=(int(n_rep), coefficients.size), replace=True)
    surrogates = (signs * coefficients[None, :]) @ mem.T
    surrogates += x.mean()
    return np.asarray(surrogates, dtype=np.float64)


def validate_surrogates(
    values: np.ndarray, surrogates: np.ndarray, basis: MoranBasis
) -> dict[str, Any]:
    x = np.asarray(values, dtype=float)
    s = np.asarray(surrogates, dtype=float)
    if s.ndim != 2 or s.shape[1] != x.size:
        raise ValueError("Surrogate array has invalid shape")
    observed_i = moran_i(x, basis.weights)
    null_i = np.array([moran_i(row, basis.weights) for row in s], dtype=float)
    observed_mean = float(x.mean())
    observed_sd = float(x.std(ddof=1))
    null_mean = s.mean(axis=1)
    null_sd = s.std(axis=1, ddof=1)
    return {
        "observed_moran_i": observed_i,
        "null_moran_i_median": float(np.nanmedian(null_i)),
        "null_moran_i_min": float(np.nanmin(null_i)),
        "null_moran_i_max": float(np.nanmax(null_i)),
        "max_abs_moran_difference": float(np.nanmax(np.abs(null_i - observed_i))),
        "observed_mean": observed_mean,
        "max_abs_mean_difference": float(np.nanmax(np.abs(null_mean - observed_mean))),
        "observed_sd": observed_sd,
        "max_abs_sd_difference": float(np.nanmax(np.abs(null_sd - observed_sd))),
        "basis_reconstruction_rmse": basis.reconstruction_rmse,
    }


def _standardize_vector(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    norm = np.sqrt(np.sum(x * x))
    if not np.isfinite(norm) or norm <= 0:
        return np.full_like(x, np.nan)
    return x / norm


def pearson_correlation(x: np.ndarray, y: np.ndarray) -> float:
    xz = _standardize_vector(np.asarray(x, dtype=float))
    yz = _standardize_vector(np.asarray(y, dtype=float))
    if not np.isfinite(xz).all() or not np.isfinite(yz).all():
        return float("nan")
    return float(xz @ yz)


def spearman_correlation(x: np.ndarray, y: np.ndarray) -> float:
    return pearson_correlation(rankdata(x, method="average"), rankdata(y, method="average"))


def rowwise_pearson(matrix: np.ndarray, vector: np.ndarray) -> np.ndarray:
    x = np.asarray(matrix, dtype=float)
    y = np.asarray(vector, dtype=float).reshape(-1)
    if x.ndim != 2 or x.shape[1] != y.size:
        raise ValueError("Matrix/vector dimensions differ")
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean()
    xnorm = np.sqrt(np.sum(xc * xc, axis=1))
    ynorm = np.sqrt(np.sum(yc * yc))
    denom = xnorm * ynorm
    out = np.full(x.shape[0], np.nan, dtype=float)
    valid = denom > 0
    out[valid] = (xc[valid] @ yc) / denom[valid]
    return out


def rowwise_spearman(matrix: np.ndarray, vector: np.ndarray) -> np.ndarray:
    ranks_x = rankdata(np.asarray(matrix, dtype=float), axis=1, method="average")
    ranks_y = rankdata(np.asarray(vector, dtype=float), method="average")
    return rowwise_pearson(ranks_x, ranks_y)


def empirical_p(observed: float, null_values: np.ndarray) -> float:
    null = np.asarray(null_values, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or null.size == 0:
        return float("nan")
    count = int(np.sum(np.abs(null) >= abs(observed)))
    return float((count + 1) / (null.size + 1))


def adjust_bh(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if pv.size == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    n = ranked.size
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    inv = np.empty(n, dtype=int)
    inv[order] = np.arange(n)
    out[valid] = adjusted[inv]
    return out


def file_sha256(path: Union[str, Path], *, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def deterministic_seed(base_seed: int, *parts: Any) -> int:
    digest = hashlib.sha256(
        "|".join([str(base_seed)] + [str(x) for x in parts]).encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:4], byteorder="little", signed=False)


BRANCHES = ("LC", "SNVTA")
MODELS = ("lme_intercept",)
# Task IDs identify the corresponding seeded analyses.
TASK_IDS = {("LC", "lme_intercept"): 1, ("SNVTA", "lme_intercept"): 3}


def scoped_task_ids():
    return TASK_IDS


SETTINGS = {
    "atlas": {
        "weight_k_neighbors": 10,
        "weight_kernel": "gaussian",
        "weight_bandwidth": "median_neighbor",
        "max_weight_components": 1,
    },
    "scope": "cortex_muse",
    "coverage_threshold": 0.70,
    "min_regions": 40,
    "n_null": 10000,
    "chunk": 500,
    "base_seed": 20260803,
    "metric": "spearman",
    "two_sided": True,
    "preserve_distribution": False,
    "fdr_group_by": ["map_id", "scope", "metric", "pet_analysis_family"],
    "global_plsc_gate": 0.05,
}
COMPARISON_TOLERANCES = {
    "salience_atol": 1e-10,
    "salience_rtol": 1e-8,
    "spearman_atol": 1e-12,
    "spearman_rtol": 1e-10,
    "probability_atol": 1e-12,
    "probability_rtol": 1e-10,
    "counts_and_order": "exact",
    "monte_carlo": "No Monte Carlo waiver with unchanged ordered maps, seed chunks and runtime; changed inputs remain conditional comparisons",
}


def stable_hash(value, *, length=64):
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            default=lambda v: v.item() if isinstance(v, np.generic) else str(v),
        ).encode()
    ).hexdigest()[:length]


def write_json(value, path):
    path = Path(path)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(
        json.dumps(
            value,
            indent=2,
            sort_keys=True,
            default=lambda v: v.item() if isinstance(v, np.generic) else str(v),
        )
        + "\n"
    )
    temporary.replace(path)


def write_table(value, path):
    path = Path(path)
    temporary = path.with_name(path.name + ".partial")
    value.to_csv(temporary, sep=",", index=False, float_format="%.17g", na_rep="NA")
    temporary.replace(path)


def checked_output(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path):
    return pd.read_csv(path, sep=",", low_memory=False, dtype={"region_id": str})


def identity(paths):
    return {str(Path(p).resolve()): file_sha256(p) for p in paths}


def verify_identity(hashes):
    for path, expected in hashes.items():
        if file_sha256(path) != expected:
            raise RuntimeError("Input/code identity changed: " + path)


def load_resources(data_root):
    root = Path(data_root)
    paths = {
        k: root / f
        for k, f in {
            "atlas": "regions.csv",
            "values": "pet_regions.csv",
            "distances": "cortical_distances.csv",
            "references": "reference_values.csv",
        }.items()
    }
    atlas, values, distance_df = [read_table(paths[k]) for k in ("atlas", "values", "distances")]
    if "region_id" not in atlas:
        atlas["region_id"] = atlas.volume_index.astype(str)
    atlas.region_id = atlas.region_id.astype(str)
    if len(atlas) != 145 or atlas.region_id.duplicated().any():
        raise ValueError("Exactly 145 ordered atlas regions required")
    metadata = [
        "pet_id",
        "target",
        "tracer",
        "source",
        "pet_analysis_family",
        "family_order",
        "eligible_branches",
        "allowed_scopes",
        "qc_status",
    ]
    if not set(metadata + ["region_id", "pet_value", "coverage"]) <= set(values):
        raise ValueError("PET resource lacks required measurement/metadata fields")
    pets = values[metadata].drop_duplicates()
    if (
        len(pets) != 22
        or pets.pet_id.duplicated().any()
        or pets.pet_analysis_family.value_counts().to_dict() != {"exploratory": 16, "primary": 6}
    ):
        raise ValueError("22 frozen PET maps with consistent metadata and 6/16 families required")
    if (
        not pets.qc_status.isin(["PASS", "PASS_WITH_FLAG"]).all()
        or values.duplicated(["pet_id", "region_id"]).any()
    ):
        raise ValueError("PET QC/identity invalid")
    ids = atlas.loc[
        atlas.in_cortex_muse.astype(str).str.lower().isin(["true", "1"]), "region_id"
    ].tolist()
    if (
        distance_df.iloc[:, 0].astype(str).tolist() != ids
        or distance_df.columns[1:].astype(str).tolist() != ids
    ):
        raise ValueError("Distance matrix row/column order must equal cortical atlas order")
    distance = distance_df.iloc[:, 1:].to_numpy(float)
    if (
        not np.isfinite(distance).all()
        or np.any(distance < 0)
        or not np.allclose(distance, distance.T, atol=1e-8, rtol=1e-7)
    ):
        raise ValueError("Invalid frozen distances")
    np.fill_diagonal(distance, 0)
    return {
        "atlas": atlas,
        "scope_ids": ids,
        "pets": pets,
        "values": values,
        "distance": distance,
        "identity": identity(paths.values()),
        "paths": paths,
    }


def seed_schedule(data_root):
    # Read the fixed scientific chunk-seed schedule.
    r = read_table(Path(data_root) / "reference_values.csv")
    r = r[r.scope.eq("molecular_seed")].copy()
    required = [
        "map_id",
        "molecular_scope",
        "perm_start",
        "perm_end",
        "seed",
        "config_hash",
        "analysis_fingerprint",
    ]
    if not set(required) <= set(r):
        raise ValueError("Original scientific seed metadata missing")
    r = r[required].rename(columns={"molecular_scope": "scope"})
    for k in ("perm_start", "perm_end", "seed"):
        r[k] = pd.to_numeric(r[k], errors="raise").astype("int64")
    return r


def runtime():
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "executable": sys.executable,
        "thread_settings": {
            k: os.getenv(k, "")
            for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]
        },
    }


def scientific_map_id(branch, model):
    return f"primary_dxbl_all__{branch}__{model}__LV1"


def work_receipt(results, task_id, value):
    work = checked_output(Path(results) / ".work" / "molecular")
    required = []
    if value.get("model_status") == "estimated":
        required = [
            str(Path(value["output_file"]).resolve()),
            str(
                (
                    Path(results) / ".work" / "molecular" / "tasks" / f"task_{task_id:03d}.json"
                ).resolve()
            ),
        ]
    write_json(
        {
            **value,
            "required_files": required,
            "signature": os.environ.get("MODEL_SIGNATURE", ""),
            "profile": os.environ.get("MODEL_PROFILE", ""),
        },
        work / f"task_{task_id:04d}.json",
    )


def dispatch_ids(config):
    selected = set(config["selected_task_ids"])
    if not selected <= set(scoped_task_ids().values()):
        raise ValueError("Prepared molecular tasks are outside this release scope")
    if os.environ.get("FIGURE4_TASK_IDS"):
        requested = {int(i) for i in os.environ["FIGURE4_TASK_IDS"].split(",")}
        if not requested or not requested <= selected:
            raise ValueError(
                "FIGURE4_TASK_IDS includes an unapproved or unprepared molecular task"
            )
        return sorted(requested)
    return sorted(selected)


def refresh_pending(results, out, config):
    work = checked_output(Path(results) / ".work" / "molecular")
    pending = []
    for task_id in dispatch_ids(config):
        status_path = work / f"task_{task_id:04d}.json"
        if not status_path.exists():
            pending.append(task_id)
            continue
        receipt = json.loads(status_path.read_text())
        if (
            receipt.get("signature") != config["launcher_signature"]
            or receipt.get("analysis_fingerprint") != config["analysis_fingerprint"]
        ):
            raise RuntimeError("Completed molecular receipt has another task/run identity")
        if receipt.get("status") != "complete":
            raise RuntimeError("Existing molecular task receipt is not complete")
        if receipt.get("output_file") and receipt.get("output_sha256") != file_sha256(
            receipt["output_file"]
        ):
            raise RuntimeError("Completed molecular null output checksum changed")
    (work / "n_tasks.txt").write_text(str(max(TASK_IDS.values())) + "\n")
    (work / "pending_task_ids.txt").write_text("".join(str(i) + "\n" for i in pending))
    write_json(
        {
            "signature": config["launcher_signature"],
            "profile": config["launcher_profile"],
            "analysis_fingerprint": config["analysis_fingerprint"],
            "selected_task_ids": config["selected_task_ids"],
            "prepared_file": str(out / "prepared.json"),
            "prepared_sha256": file_sha256(out / "prepared.json"),
        },
        work / "preparation_receipt.json",
    )


def prepare(results, plsc_root, pet_root, *, profile="study", maps=None, n_null=None):
    out = checked_output(Path(results) / ".work" / "molecular")
    if (out / "prepared.json").exists():
        out, previous = load_prepared(results)
        requested = (
            sorted(maps) if maps else sorted(scientific_map_id(b, m) for b, m in scoped_task_ids())
        )
        if (
            previous["requested_maps"] != requested
            or previous["profile"] != profile
            or previous["settings"]["n_null"] != (n_null or 10000)
        ):
            raise RuntimeError("Molecular resume selection/profile/draw count changed")
        refresh_pending(results, out, previous)
        return {
            "status": "prepared_unchanged",
            "n_tasks": previous["n_tasks"],
            "n_observed": previous["n_observed"],
        }
    resources = load_resources(pet_root)
    seed_file = Path(pet_root) / "reference_values.csv"
    historical_seeds = seed_schedule(pet_root)
    seed_keys = ["map_id", "scope", "perm_start", "perm_end"]
    if historical_seeds.duplicated(seed_keys).any():
        raise ValueError("Historical scientific seed manifest has duplicate chunk identities")
    settings = dict(SETTINGS)
    if profile == "study" and n_null not in (None, 10000):
        raise ValueError("Study molecular results require 10000 nulls")
    settings["n_null"] = 10000 if n_null is None else int(n_null)
    if settings["n_null"] < 1:
        raise ValueError("At least one null draw required")
    global_file = Path(plsc_root) / "plsc_global_results.csv"
    region_file = Path(plsc_root) / "plsc_regional_results.csv"
    global_results = read_table(global_file).rename(columns={"P": "permutation_p"})
    global_results["status"] = np.where(
        np.isfinite(global_results.permutation_p), "estimated", "unestimated"
    )
    region_results = read_table(region_file)
    required = {"branch", "model", "lv", "permutation_p", "status"}
    if not required <= set(global_results) or not {
        "branch",
        "model",
        "lv",
        "region_id",
        "salience_oriented",
    } <= set(region_results):
        raise ValueError("PLSC molecular interface is incomplete")
    selected = global_results[
        global_results.branch.isin(BRANCHES)
        & global_results.model.isin(MODELS)
        & (global_results.lv == "LV1")
    ].copy()
    selected["map_id"] = [scientific_map_id(b, m) for b, m in zip(selected.branch, selected.model)]
    if selected.map_id.duplicated().any():
        raise ValueError("Duplicate PLSC map scientific identity")
    requested = set(maps) if maps else {scientific_map_id(b, m) for b, m in scoped_task_ids()}
    if not requested <= {scientific_map_id(b, m) for b, m in scoped_task_ids()}:
        raise ValueError(
            "Unapproved molecular map requested; slope, common-X and JOINT_RAW2D are excluded"
        )
    missing = requested - set(selected.map_id)
    if missing:
        raise ValueError(
            "Required current PLSC map records are absent: " + ",".join(sorted(missing))
        )
    selected = selected[selected.map_id.isin(requested)].sort_values("map_id")
    if profile == "study":
        if not {"resampling_profile", "num_permutations", "num_bootstraps"} <= set(selected):
            raise ValueError(
                "Study molecular preparation requires explicit completed PLSC resampling metadata"
            )
        if not (
            selected.resampling_profile.eq("study")
            & selected.num_permutations.eq(5000)
            & selected.num_bootstraps.eq(5000)
        ).all():
            raise ValueError(
                "Developmental PLSC results cannot enter final study molecular inference"
            )
    source_identity = {
        **resources["identity"],
        **identity([global_file, region_file, seed_file, __file__]),
    }
    config_hash = stable_hash(settings, length=24)
    fingerprint = stable_hash(
        {
            "input_code_identity": source_identity,
            "profile": profile,
            "requested": sorted(requested),
            "settings": settings,
        }
    )
    observations, exclusions, coordinate_rows, task_rows = [], [], [], []
    pet_lookup = {p: g.set_index("region_id") for p, g in resources["values"].groupby("pet_id")}
    atlas_index = resources["atlas"].set_index("region_id")
    (out / "matrices").mkdir()
    (out / "tasks").mkdir()
    for row in selected.to_dict("records"):
        map_id, branch, model = row["map_id"], row["branch"], row["model"]
        task_id = TASK_IDS[(branch, model)]
        p = float(row["permutation_p"])
        if np.isfinite(p) and not 0 <= p <= 1:
            raise ValueError("Invalid global PLSC probability: " + map_id)
        if (
            row["status"] not in ("estimated", "complete", "COMPLETE", "ok", "OK")
            or not np.isfinite(p)
            or p >= 0.05
        ):
            exclusions.append(
                {
                    "map_id": map_id,
                    "pet_id": "",
                    "scope": settings["scope"],
                    "reason": "plsc_not_estimated_or_global_permutation_p_not_below_0.05",
                }
            )
            continue
        current = region_results[
            (region_results.branch == branch)
            & (region_results.model == model)
            & (region_results.lv == "LV1")
        ]
        if (
            len(current) != 145
            or current.region_id.duplicated().any()
            or set(current.region_id) != set(resources["atlas"].region_id)
        ):
            raise ValueError(
                "PLSC map must retain all 145 distinct atlas region identities: " + map_id
            )
        series = current.set_index("region_id").salience_oriented.astype(float)
        ids = [rid for rid in resources["scope_ids"] if np.isfinite(series.loc[rid])]
        if len(ids) < settings["min_regions"]:
            exclusions.append(
                {
                    "map_id": map_id,
                    "pet_id": "",
                    "scope": settings["scope"],
                    "reason": "map_has_too_few_finite_regions",
                }
            )
            continue
        x = series.loc[ids].to_numpy(float)
        positions = [resources["scope_ids"].index(r) for r in ids]
        weights, weight_meta = build_spatial_weights(
            resources["distance"][np.ix_(positions, positions)], settings
        )
        basis = compute_moran_basis(weights)
        matrix_path = out / "matrices" / (map_id + ".npz")
        np.savez_compressed(
            matrix_path,
            region_ids=np.asarray(ids),
            map_values=x,
            weights=basis.weights,
            eigenvalues=basis.eigenvalues,
            eigenvectors=basis.eigenvectors,
            reconstruction_rmse=np.asarray(basis.reconstruction_rmse),
        )
        before = len(observations)
        for pet in resources["pets"].to_dict("records"):
            pet_id = pet["pet_id"]
            reason = ""
            if branch not in str(pet["eligible_branches"]).split(","):
                reason = "pet_not_eligible_for_branch"
            elif settings["scope"] not in str(pet["allowed_scopes"]).split(","):
                reason = "pet_allowed_scopes_not_compatible"
            pet_values = pet_lookup[pet_id]
            pair_ids = [
                r
                for r in ids
                if r in pet_values.index
                and np.isfinite(pet_values.loc[r, "pet_value"])
                and np.isfinite(pet_values.loc[r, "coverage"])
                and pet_values.loc[r, "coverage"] >= settings["coverage_threshold"]
            ]
            if len(pair_ids) < settings["min_regions"]:
                reason = (
                    reason or f'pair_has_too_few_regions:{len(pair_ids)}<{settings["min_regions"]}'
                )
            xx, yy = series.loc[pair_ids].to_numpy(float), pet_values.loc[
                pair_ids, "pet_value"
            ].to_numpy(float)
            if len(pair_ids) and (np.std(xx, ddof=1) <= 0 or np.std(yy, ddof=1) <= 0):
                reason = reason or "constant_map"
            if reason:
                exclusions.append(
                    {
                        "map_id": map_id,
                        "pet_id": pet_id,
                        "scope": settings["scope"],
                        "reason": reason,
                    }
                )
                continue
            pair_id = f'{map_id}__{settings["scope"]}__{pet_id}'
            analysis_id = pair_id + "__spearman"
            observations.append(
                {
                    **pet,
                    "analysis_id": analysis_id,
                    "pair_id": pair_id,
                    "map_id": map_id,
                    "branch": branch,
                    "model": model,
                    "lv": "LV1",
                    "scope": settings["scope"],
                    "metric": "spearman",
                    "map_statistic": "salience",
                    "lv_permutation_p": p,
                    "effect": spearman_correlation(xx, yy),
                    "n_regions": len(pair_ids),
                    "mask_hash": stable_hash(pair_ids, length=20),
                    "mask_policy": "pairwise_fixed",
                    "fdr_family_id": f'{map_id}__{settings["scope"]}__spearman__{pet["pet_analysis_family"]}',
                }
            )
            coordinate_rows.extend(
                {
                    "analysis_id": analysis_id,
                    "map_id": map_id,
                    "pet_id": pet_id,
                    "region_id": rid,
                    "position": i,
                    "map_position": ids.index(rid),
                    "region_name": atlas_index.loc[rid, "region_name"],
                    "hemisphere": atlas_index.loc[rid, "hemisphere"],
                    "salience": float(xx[i]),
                    "pet_value": float(yy[i]),
                    "coverage": float(pet_values.loc[rid, "coverage"]),
                }
                for i, rid in enumerate(pair_ids)
            )
        if len(observations) > before:
            task_rows.append(
                {
                    "task_id": task_id,
                    "map_id": map_id,
                    "matrix_file": str(matrix_path),
                    "n_null": settings["n_null"],
                    "analysis_fingerprint": fingerprint,
                    **weight_meta,
                }
            )
    columns = ["map_id", "pet_id", "scope", "reason"]
    observed = pd.DataFrame(observations)
    if observed.empty:
        observed = pd.DataFrame(
            columns=[
                "map_id",
                "branch",
                "model",
                "pet_id",
                "metric",
                "scope",
                "analysis_id",
                "pet_analysis_family",
                "effect",
            ]
        )
    if not observed.empty:
        observed = observed.sort_values(["map_id", "metric", "pet_id", "analysis_id"]).reset_index(
            drop=True
        )
    write_table(observed, out / "molecular_observed.csv")
    write_table(pd.DataFrame(exclusions, columns=columns), out / "molecular_unestimated.csv")
    write_table(
        pd.DataFrame(
            coordinate_rows,
            columns=[
                "analysis_id",
                "map_id",
                "pet_id",
                "region_id",
                "position",
                "map_position",
                "region_name",
                "hemisphere",
                "salience",
                "pet_value",
                "coverage",
            ],
        ),
        out / "molecular_plot_coordinates.csv",
    )
    write_table(
        pd.DataFrame(
            task_rows,
            columns=[
                "task_id",
                "map_id",
                "matrix_file",
                "n_null",
                "analysis_fingerprint",
                "weight_kernel",
                "weight_k_neighbors",
                "weight_bandwidth",
                "weight_n_components",
                "weight_sum",
                "weight_hash",
            ],
        ),
        out / "molecular_tasks.csv",
    )
    frozen_outputs = [
        out / x
        for x in [
            "molecular_observed.csv",
            "molecular_unestimated.csv",
            "molecular_plot_coordinates.csv",
            "molecular_tasks.csv",
        ]
    ]
    frozen_outputs += list((out / "matrices").glob("*.npz"))
    chunks = []
    for task in task_rows:
        for start in range(0, settings["n_null"], settings["chunk"]):
            end = min(settings["n_null"], start + settings["chunk"])
            # Scientific seeds are independent of output paths and task batching.
            old = historical_seeds[
                (historical_seeds.map_id == task["map_id"])
                & (historical_seeds.scope == settings["scope"])
                & (historical_seeds.perm_start == start)
            ]
            if len(old) != 1 or int(old.iloc[0].perm_end) < end:
                raise ValueError("Required original map/chunk scientific seed is unavailable")
            old = old.iloc[0]
            expected_seed = deterministic_seed(
                settings["base_seed"],
                task["map_id"],
                settings["scope"],
                start,
                int(old.perm_end),
                old.config_hash,
                old.analysis_fingerprint,
            )
            if int(old.seed) != expected_seed:
                raise ValueError(
                    "Recorded original seed does not reproduce the original hash procedure"
                )
            chunks.append(
                {
                    "task_id": task["task_id"],
                    "map_id": task["map_id"],
                    "scope": settings["scope"],
                    "perm_start": start,
                    "perm_end": end,
                    "seed": int(old.seed),
                    "seed_reference_config_hash": old.config_hash,
                    "seed_reference_input_fingerprint": old.analysis_fingerprint,
                    "seed_reference_perm_end": int(old.perm_end),
                }
            )
    write_table(
        pd.DataFrame(
            chunks,
            columns=[
                "task_id",
                "map_id",
                "scope",
                "perm_start",
                "perm_end",
                "seed",
                "seed_reference_config_hash",
                "seed_reference_input_fingerprint",
                "seed_reference_perm_end",
            ],
        ),
        out / "molecular_null_chunks.csv",
    )
    frozen_outputs.append(out / "molecular_null_chunks.csv")
    config = {
        "status": "prepared",
        "execution_scope": "main",
        "profile": profile,
        "settings": settings,
        "config_hash": config_hash,
        "analysis_fingerprint": fingerprint,
        "input_code_identity": source_identity,
        "prepared_identity": identity(frozen_outputs),
        "runtime": runtime(),
        "requested_maps": sorted(requested),
        "n_tasks": len(task_rows),
        "n_observed": len(observed),
        "selected_task_ids": sorted(
            TASK_IDS[(r["branch"], r["model"])] for r in selected.to_dict("records")
        ),
        "launcher_signature": os.environ.get("MODEL_SIGNATURE", ""),
        "launcher_profile": os.environ.get("MODEL_PROFILE", ""),
        "comparison_tolerances": COMPARISON_TOLERANCES,
        "plsc_regional_file": str(region_file.resolve()),
        "data_root": str(Path(pet_root).resolve()),
        "bootstrap_available": False,
        "bootstrap_reason": "Original molecular scope has no included regional bootstrap replicate analysis",
        "seed_policy": "Recorded original scientific map/chunk seeds, revalidated by original deterministic_seed; new independent input/code/config identity governs resume",
    }
    write_json(config, out / "prepared.json")
    estimable_ids = {r["task_id"] for r in task_rows}
    for r in selected.to_dict("records"):
        task_id = TASK_IDS[(r["branch"], r["model"])]
        if task_id not in estimable_ids:
            work_receipt(
                results,
                task_id,
                {
                    "status": "complete",
                    "model_status": "not_estimable",
                    "task_id": task_id,
                    "map_id": r["map_id"],
                    "analysis_fingerprint": fingerprint,
                    "reason": ";".join(
                        sorted({x["reason"] for x in exclusions if x["map_id"] == r["map_id"]})
                    ),
                },
            )
    refresh_pending(results, out, config)
    return {
        "status": "prepared",
        "n_tasks": len(task_rows),
        "n_observed": len(observed),
        "output": str(out),
    }


def load_prepared(results):
    out = checked_output(Path(results) / ".work" / "molecular")
    config = json.loads((out / "prepared.json").read_text())
    if config["launcher_signature"] != os.environ.get("MODEL_SIGNATURE", "") or config[
        "launcher_profile"
    ] != os.environ.get("MODEL_PROFILE", ""):
        raise RuntimeError("Molecular launcher signature/profile changed")
    if config.get("execution_scope", "full") != "main":
        raise RuntimeError("Molecular execution scope changed")
    verify_identity(config["input_code_identity"])
    verify_identity(config["prepared_identity"])
    if runtime() != config["runtime"]:
        raise RuntimeError("Molecular Python runtime or thread settings changed")
    return out, config


def task(results, task_id):
    out, config = load_prepared(results)
    tasks = read_table(out / "molecular_tasks.csv")
    selected = tasks[tasks.task_id == task_id]
    if len(selected) != 1:
        receipt_path = out / f"task_{task_id:04d}.json"
        if receipt_path.is_file():
            r = json.loads(receipt_path.read_text())
            if (
                r.get("status") == "complete"
                and r.get("model_status") == "not_estimable"
                and r.get("signature") == config["launcher_signature"]
                and r.get("analysis_fingerprint") == config["analysis_fingerprint"]
            ):
                return {
                    "status": "complete",
                    "model_status": "not_estimable",
                    "task_id": task_id,
                    "reason": r["reason"],
                }
        raise ValueError(
            "Task ID does not identify one molecular map or explicit unestimated receipt"
        )
    definition = selected.iloc[0].to_dict()
    status_path = out / "tasks" / f"task_{task_id:03d}.json"
    npz_path = out / "tasks" / f"task_{task_id:03d}.npz"
    if status_path.exists():
        status = json.loads(status_path.read_text())
        if (
            status.get("status") == "complete"
            and status.get("analysis_fingerprint") == config["analysis_fingerprint"]
            and status.get("output_sha256") == file_sha256(npz_path)
        ):
            work_receipt(
                results,
                task_id,
                {
                    "status": "complete",
                    "model_status": "estimated",
                    "task_id": task_id,
                    "map_id": definition["map_id"],
                    "analysis_fingerprint": config["analysis_fingerprint"],
                    "output_file": str(npz_path),
                    "output_sha256": file_sha256(npz_path),
                },
            )
            return {"status": "resumed_unchanged", "task_id": task_id}
        raise RuntimeError("Existing molecular result has an invalid identity")
    start_time = time.monotonic()
    with np.load(definition["matrix_file"], allow_pickle=False) as matrix:
        ids, x = matrix["region_ids"].astype(str), matrix["map_values"]
        basis = MoranBasis(
            matrix["eigenvectors"],
            matrix["eigenvalues"],
            matrix["weights"],
            float(matrix["reconstruction_rmse"]),
        )
    observed = read_table(out / "molecular_observed.csv")
    observed = observed[observed.map_id == definition["map_id"]].sort_values(
        ["metric", "pet_id", "analysis_id"]
    )
    coordinates = read_table(out / "molecular_plot_coordinates.csv")
    coordinate_lookup = {
        aid: group.sort_values("position") for aid, group in coordinates.groupby("analysis_id")
    }
    chunks = read_table(out / "molecular_null_chunks.csv")
    chunks = chunks[chunks.task_id == task_id].sort_values("perm_start")
    null_matrix = np.empty((config["settings"]["n_null"], len(observed)))
    qcs = []
    for chunk in chunks.to_dict("records"):
        begin, end = int(chunk["perm_start"]), int(chunk["perm_end"])
        surrogates = moran_randomize(
            x,
            basis,
            n_rep=end - begin,
            random_state=int(chunk["seed"]),
        )
        qc = validate_surrogates(x, surrogates, basis)
        qc["pass"] = (
            qc["max_abs_moran_difference"] <= 1e-7
            and qc["max_abs_mean_difference"] <= 1e-9
            and qc["max_abs_sd_difference"] <= 1e-8
            and qc["basis_reconstruction_rmse"] <= 1e-9
        )
        qcs.append({**chunk, **qc})
        if not qc["pass"]:
            write_json({"status": "failed_null_qc", "task_id": task_id, "qc": qcs}, status_path)
            raise RuntimeError("Moran surrogate QC failed; no inference finalized")
        for column, obs in enumerate(observed.to_dict("records")):
            coords = coordinate_lookup[obs["analysis_id"]]
            positions = coords.map_position.to_numpy(int)
            if not np.array_equal(ids[positions], coords.region_id.to_numpy(str)):
                raise RuntimeError("Molecular pair mask and null parcel ordering differ")
            null_matrix[begin:end, column] = rowwise_spearman(
                surrogates[:, positions], coords.pet_value.to_numpy(float)
            )
    if not np.isfinite(null_matrix).all():
        raise RuntimeError("Null correlations contain nonfinite values")
    np.savez_compressed(
        npz_path,
        null_values=null_matrix,
        analysis_ids=observed.analysis_id.to_numpy(str),
        region_ids=ids,
        analysis_fingerprint=np.asarray(config["analysis_fingerprint"]),
    )
    write_json(
        {
            "status": "complete",
            "task_id": task_id,
            "map_id": definition["map_id"],
            "analysis_fingerprint": config["analysis_fingerprint"],
            "runtime": runtime(),
            "host": platform.node(),
            "job_id": os.getenv("JOB_ID", ""),
            "sge_task_id": os.getenv("SGE_TASK_ID", ""),
            "wall_seconds": time.monotonic() - start_time,
            "output_sha256": file_sha256(npz_path),
            "n_null": len(null_matrix),
            "n_analyses": len(observed),
            "chunk_qc": qcs,
        },
        status_path,
    )
    work_receipt(
        results,
        task_id,
        {
            "status": "complete",
            "model_status": "estimated",
            "task_id": task_id,
            "map_id": definition["map_id"],
            "analysis_fingerprint": config["analysis_fingerprint"],
            "output_file": str(npz_path),
            "output_sha256": file_sha256(npz_path),
        },
    )
    return {
        "status": "complete",
        "task_id": task_id,
        "n_null": len(null_matrix),
        "n_analyses": len(observed),
    }


def finalize(results_root):
    results = results_root
    out, config = load_prepared(results)
    refresh_pending(results, out, config)
    if (Path(results) / ".work/molecular/pending_task_ids.txt").read_text().strip():
        raise RuntimeError(
            "Selected molecular tasks are incomplete; no partial-family finalization"
        )
    tasks = read_table(out / "molecular_tasks.csv")
    selected_ids = dispatch_ids(config)
    tasks = tasks[tasks.task_id.isin(selected_ids)]
    complete_scope = set(selected_ids) == set(scoped_task_ids().values())
    observed = read_table(out / "molecular_observed.csv")
    rows, qcs = [], []
    for task_row in tasks.to_dict("records"):
        task_id = int(task_row["task_id"])
        status = json.loads((out / "tasks" / f"task_{task_id:03d}.json").read_text())
        npz = out / "tasks" / f"task_{task_id:03d}.npz"
        if (
            status.get("status") != "complete"
            or status["analysis_fingerprint"] != config["analysis_fingerprint"]
            or status["output_sha256"] != file_sha256(npz)
        ):
            raise RuntimeError("Molecular task missing, incomplete, or identity mismatch")
        qcs.extend(status["chunk_qc"])
        with np.load(npz, allow_pickle=False) as saved:
            nulls, aids = saved["null_values"], saved["analysis_ids"].astype(str)
        obs = observed[observed.map_id == task_row["map_id"]].sort_values(
            ["metric", "pet_id", "analysis_id"]
        )
        if (
            not np.array_equal(aids, obs.analysis_id.to_numpy(str))
            or nulls.shape != (config["settings"]["n_null"], len(obs))
            or not np.isfinite(nulls).all()
        ):
            raise RuntimeError("Completed molecular null array coverage or order mismatch")
        if not all(q["pass"] for q in status["chunk_qc"]):
            raise RuntimeError("Unresolved null QC failure")
        for column, row in enumerate(obs.to_dict("records")):
            rows.append(
                {
                    **row,
                    "p_spatial": empirical_p(row["effect"], nulls[:, column]),
                    "n_null": len(nulls),
                    "status": "estimated",
                    "profile": config["profile"],
                }
            )
    results = pd.DataFrame(
        rows, columns=list(observed.columns) + ["p_spatial", "n_null", "status", "profile"]
    )
    results["q_bh_within_plsc_map_family"] = np.nan
    results["fdr_family_size"] = 0
    for _, index in results.groupby(SETTINGS["fdr_group_by"]).groups.items():
        expected = 6 if results.loc[index, "pet_analysis_family"].iloc[0] == "primary" else 16
        if len(index) == expected:
            results.loc[index, "q_bh_within_plsc_map_family"] = adjust_bh(
                results.loc[index, "p_spatial"].to_numpy(float)
            )
        results.loc[index, "fdr_family_size"] = len(index)
    results["final_inference"] = (
        config["profile"] == "study" and config["settings"]["n_null"] == 10000
    )
    results["complete_approved_map_coverage"] = complete_scope
    # Keep all 44 planned map pairs. Missing values remain NA; never shrink BH.
    planned = []
    resources = load_resources(config["data_root"])
    for branch in BRANCHES:
        map_id = scientific_map_id(branch, "lme_intercept")
        for pet in resources["pets"].to_dict("records"):
            found = results[(results.map_id == map_id) & (results.pet_id == pet["pet_id"])]
            planned.append(
                found.iloc[0].to_dict()
                if len(found)
                else {
                    **pet,
                    "map_id": map_id,
                    "branch": branch,
                    "model": "lme_intercept",
                    "effect": np.nan,
                    "p_spatial": np.nan,
                    "q_bh_within_plsc_map_family": np.nan,
                }
            )
    publication = pd.DataFrame(planned)
    # A family is final only if every planned test was estimable and accounted for.
    for (_, family), idx in publication.groupby(["map_id", "pet_analysis_family"]).groups.items():
        expected = 6 if family == "primary" else 16
        if len(idx) != expected:
            raise RuntimeError("Planned molecular family changed")
        if not np.isfinite(publication.loc[idx, "p_spatial"]).all():
            publication.loc[idx, "q_bh_within_plsc_map_family"] = np.nan
    publication = publication.rename(
        columns={"effect": "rho", "p_spatial": "P", "q_bh_within_plsc_map_family": "P_fdr"}
    )
    keep = [
        "branch",
        "model",
        "map_id",
        "pet_id",
        "target",
        "tracer",
        "source",
        "pet_analysis_family",
        "family_order",
        "rho",
        "P",
        "P_fdr",
        "n_regions",
        "n_null",
    ]
    dest = checked_output(Path(results_root) / "figure4")
    write_table(
        publication.reindex(columns=keep).sort_values(
            ["map_id", "pet_analysis_family", "family_order", "pet_id"]
        ),
        dest / "molecular_results.csv",
    )
    # A gated parent has no null draws; retain a declared empty diagnostic schema.
    empty_qc_columns = [
        "basis_reconstruction_rmse",
        "map_id",
        "max_abs_mean_difference",
        "max_abs_moran_difference",
        "max_abs_sd_difference",
        "null_moran_i_max",
        "null_moran_i_median",
        "null_moran_i_min",
        "observed_mean",
        "observed_moran_i",
        "observed_sd",
        "pass",
        "perm_end",
        "perm_start",
        "scope",
        "seed",
        "seed_reference_config_hash",
        "seed_reference_input_fingerprint",
        "seed_reference_perm_end",
        "task_id",
    ]
    write_table(
        pd.DataFrame(qcs) if qcs else pd.DataFrame(columns=empty_qc_columns),
        out / "null_diagnostics.csv",
    )
    write_table(read_table(out / "molecular_unestimated.csv"), out / "diagnostics.csv")
    return {
        "status": "complete",
        "planned_pairs": 44,
        "estimated_pairs": int(np.isfinite(publication.P).sum()),
    }


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--mode", choices=["prepare", "task", "finalize"], required=True)
    a.add_argument("--data", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--task-id", type=int)
    a.add_argument("--synthetic-n-null", type=int)
    args = a.parse_args()
    if args.synthetic_n_null is not None:
        if not (Path(args.data) / "SYNTHETIC.txt").is_file():
            raise ValueError("Reduced draws require invented fixture marker")
        os.environ["MODEL_PROFILE"] = "selftest"
    if args.mode == "prepare":
        result = prepare(
            args.out,
            Path(args.out) / "figure4",
            args.data,
            profile="selftest" if args.synthetic_n_null else "study",
            n_null=args.synthetic_n_null,
        )
    elif args.mode == "task":
        result = task(args.out, args.task_id)
    else:
        result = finalize(args.out)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
