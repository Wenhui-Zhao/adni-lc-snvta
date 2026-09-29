#!/usr/bin/env bash
# One numerical launcher. The embedded coordinator uses Python's standard library.
set -euo pipefail
REPLICATION_CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${REPLICATION_SETTINGS:-}" ]]; then source "$REPLICATION_SETTINGS"; fi
export REPLICATION_CODE_ROOT
export SGE_COMPLETION_WAIT_SECONDS="${SGE_COMPLETION_WAIT_SECONDS:-300}"
if [[ "${1:-}" == export ]]; then
  shift
  exec "${PYTHON_BIN:-python3}" "$REPLICATION_CODE_ROOT/table_views.py" "$@"
fi
exec "${PYTHON_BIN:-python3}" - "$@" <<'PY'
import argparse, csv, hashlib, json, os, pathlib, shlex, shutil, subprocess, sys, time

P = pathlib.Path
code = P(os.environ["REPLICATION_CODE_ROOT"])
a = argparse.ArgumentParser(
    description="Main-only, 14 prepared inputs to 11 public CSV views; no installations."
)
a.add_argument("command", choices=["check", "all", "figure1", "figure2", "figure3",
                                   "figure4", "figure5", "stage", "_worker"])
a.add_argument("--data", required=True)
a.add_argument("--out", required=True)
a.add_argument("--mode", choices=["local", "sge"], default="local")
a.add_argument("--stage", default="all")
a.add_argument("--task", type=int)
a.add_argument("--jobs", type=int, default=2)
x = a.parse_args()
if not 1 <= x.jobs <= 128:
    a.error("--jobs must be 1..128; local is always sequential")
try:
    completion_wait = float(os.environ.get("SGE_COMPLETION_WAIT_SECONDS", "300"))
except ValueError:
    a.error("SGE_COMPLETION_WAIT_SECONDS must be a finite positive number")
if not 0 < completion_wait < float("inf"):
    a.error("SGE_COMPLETION_WAIT_SECONDS must be a finite positive number")
data = P(x.data).resolve(strict=True)
run_root = P(x.out).resolve()
out = run_root / "private" / "internal"
if not data.is_dir():
    a.error("--data must be a directory")
if run_root == data or data in run_root.parents or run_root in data.parents:
    a.error("Input and output roots must be separate")
if run_root == code or code in run_root.parents:
    a.error("Generated outputs must be outside source")
out.mkdir(parents=True, exist_ok=True)
work = out / ".work"
work.mkdir(exist_ok=True)
for name in ["tmp", "logs", "receipts"]:
    (work / name).mkdir(exist_ok=True)
os.environ.update(
    ADNI_MODEL_DATA_DIR=str(data),
    ADNI_MODEL_RESULTS_DIR=str(out),
    MODEL_PROFILE="analysis",
    MODEL_EXECUTION_SCOPE="main",
    TMPDIR=str(work / "tmp"),
    TMP=str(work / "tmp"),
    TEMP=str(work / "tmp"),
    PYTHONDONTWRITEBYTECODE="1",
    XDG_CACHE_HOME=str(work / "cache"),
    MATLAB_PREFDIR=str(work / "matlab_prefs"),
    MCR_CACHE_ROOT=str(work / "mcr_cache"),
)
for name in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
]:
    os.environ[name] = "1"
# Production counts are fixed; synthetic tests invoke numerical functions separately.
for name in ["R_HOME", "RHOME"]:
    os.environ.pop(name, None)
rscript = os.environ.get("RSCRIPT_BIN", "Rscript")
matlab = os.environ.get("MATLAB_BIN", "matlab")
toolbox = os.environ.get("PLS_TOOLBOX_ROOT", "")
def stage(group, script, inputs, tasks, packages, memory, hours, *,
          deps=(), flags=(), final=(), prepare=(), engine="R"):
    """Execution and output ownership for one calculation stage."""
    return dict(group=group, script=script, inputs=inputs, tasks=tasks,
                packages=packages, resources=(memory, hours), deps=deps,
                flags=flags, final=final, prepare=prepare, engine=engine)


STAGES = {
    'technical': stage(
        'figure1', 'figure1', ['scans.csv'],
        tuple(range(1, 9)), ['data.table', 'lme4'], '6G', '04:00:00',
        flags=['--unit', 'technical'],
        final=[
            'private/shared/technical_cr.csv',
            'shared/technical_coefficients.csv',
            'shared/technical_references.csv',
            '.work/diagnostics/technical.csv',
        ],
    ),
    'icc': stage(
        'figure1', 'figure1', ['scans.csv'],
        tuple(range(1, 9)), ['data.table', 'lme4'], '6G', '12:00:00',
        flags=['--unit', 'icc'],
        final=['figure1/icc.csv', '.work/diagnostics/icc.csv'],
    ),
    'pathology': stage(
        'figure1', 'figure1', ['scans.csv', 'gross_pathology.csv'],
        tuple(range(1, 5)), ['data.table', 'lme4'], '6G', '08:00:00',
        flags=['--unit', 'pathology'],
        final=['figure1/gross_coefficients.csv', '.work/diagnostics/pathology.csv'],
    ),
    'figure2': stage(
        'figure2', 'figure2', ['scans.csv'],
        tuple(range(1, 3)), ['data.table', 'lme4', 'lmerTest'], '12G', '08:00:00',
        deps=['technical'],
        final=['figure2/coefficients.csv', 'figure2/reference_values.csv', '.work/diagnostics/figure2.csv'],
    ),
    'concurrent': stage(
        'figure3', 'figure3_associations', ['outcomes.csv', 'concurrent.csv'],
        tuple(range(1, 77)), ['data.table', 'lme4', 'lmerTest'], '12G', '08:00:00',
        flags=['--unit', 'concurrent'],
        final=[
            'figure3/concurrent_associations.csv',
            'figure3/concurrent_coefficients.csv',
            '.work/diagnostics/figure3_concurrent.csv',
        ],
    ),
    'prospective': stage(
        'figure3', 'figure3_associations', ['outcomes.csv', 'prospective.csv'],
        tuple(range(1, 77)), ['data.table', 'lme4', 'lmerTest'], '12G', '08:00:00',
        flags=['--unit', 'prospective'],
        final=[
            'figure3/prospective_associations.csv',
            'figure3/prospective_coefficients.csv',
            'figure3/display_metadata.csv',
            '.work/diagnostics/figure3_prospective.csv',
        ],
    ),
    'lag': stage(
        'figure3', 'figure3_lag', ['outcomes.csv', 'lag.csv', 'concurrent.csv'],
        tuple(range(1, 77)), ['data.table', 'mgcv', 'MASS'], '24G', '12:00:00',
        final=[
            'figure3/lag_models.csv',
            'figure3/lag_curves.csv',
            'figure3/lag_intervals.csv',
            'figure3/lag_interval_summary.csv',
            'figure3/lag_support.csv',
            'figure3/lag_coefficients.csv',
            '.work/diagnostics/figure3_lag.csv',
        ],
    ),
    'trajectories': stage(
        'figure4', 'figure4_trajectories', ['scans.csv', 'muse.csv', 'regions.csv', 'reference_values.csv'],
        tuple(range(1, 295)), ['data.table', 'lme4', 'digest'], '8G', '12:00:00',
        flags=['--stage', 'trajectories'],
        final=[f"private/figure4/{name}.csv" for name in
               ("trajectory_parameters", "matrix_X", "matrix_Y", "matrix_rows", "bilateral_parameters")]
        + [f"figure4/{name}.csv" for name in
           ("trajectory_coefficients", "trajectory_variance_components", "matrix_scaling", "preparation_scaling")]
        + [".work/diagnostics/trajectories.csv", ".work/plsc_inputs/regions.csv"]
        + [f".work/plsc_inputs/{region}/lme_{parameter}/{name}.csv"
           for region in ("LC", "SNVTA") for parameter in ("intercept", "slope")
           for name in ("X", "Y", "scaling", "rows")],
        prepare=['.work/trajectories/input_checks.csv'],
    ),
    'tpv': stage(
        None, 'figure4_trajectories', ['muse.csv'],
        tuple(range(1, 3)), ['data.table', 'lme4', 'digest'], '18G', '12:00:00',
        deps=['trajectories'],
        flags=['--stage', 'tpv'],
        final=[
            'private/figure4/tpv_parameters.csv',
            'figure4/tpv_scales.csv',
            'figure4/tpv_coefficients.csv',
            'figure4/tpv_variance_components.csv',
            '.work/diagnostics/tpv.csv',
        ],
    ),
    'plsc': stage(
        'figure4', 'figure4_plsc', ['regions.csv'],
        tuple(range(1, 5)), [], '16G', '12:00:00',
        deps=['trajectories'],
        engine='MATLAB',
        final=[
            'figure4/plsc_global_results.csv',
            'figure4/plsc_regional_results.csv',
            'figure4/plsc_behavior_results.csv',
            'figure4/plsc_transforms.csv',
            'private/figure4/plsc_participant_scores.csv',
            '.work/plsc/diagnostics.csv',
            '.work/plsc/completion.json',
        ],
    ),
    'molecular': stage(
        'figure4', 'figure4_molecular', ['regions.csv', 'reference_values.csv', 'pet_regions.csv', 'cortical_distances.csv'],
        (1, 3), [], '8G', '12:00:00',
        deps=['plsc'],
        engine='Python',
        final=[
            'figure4/molecular_results.csv',
            '.work/molecular/diagnostics.csv',
            '.work/molecular/null_diagnostics.csv',
        ],
    ),
    'progression': stage(
        'figure5', 'figure5_progression', ['progression.csv', 'diagnoses.csv'],
        tuple(range(1, 5)), ['data.table', 'lme4', 'digest'], '18G', '12:00:00',
        deps=['technical'],
        final=[
            'figure5/progression_results.csv',
            'figure5/progression_coefficients.csv',
            'private/figure5/progression_model_frames.csv',
            '.work/diagnostics/progression.csv',
        ],
        prepare=[
            'private/figure5/primary_landmarks.csv',
            'private/figure5/selected_index_audit.csv',
            'figure5/progression_covariate_scales.csv',
            'shared/technical_bilateral_scales.csv',
            '.work/diagnostics/progression_preparation.csv',
        ],
    ),
    'mediation': stage(
        'figure5', 'figure5_mediation', [],
        tuple(range(1, 5)), ['data.table', 'digest', 'jsonlite', 'lavaan'], '12G', '12:00:00',
        deps=['progression', 'parameters'],
        final=[
            'figure5/mediation_effects.csv',
            'figure5/mediation_parameters.csv',
            'figure5/mediation_fit_measures.csv',
            '.work/diagnostics/mediation.csv',
        ],
        prepare=['private/figure5/mediation_model_frames.csv', 'figure5/mediation_scales.csv'],
    ),
    'amyloid': stage(
        'figure5', 'figure5_memory', ['amyloid_memory.csv'],
        tuple(range(1, 5)), ['data.table', 'lme4', 'lmerTest', 'digest', 'jsonlite'], '12G', '12:00:00',
        deps=['parameters'],
        final=[
            'figure5/memory_coefficients.csv',
            'figure5/memory_variance_components.csv',
            'figure5/memory_focal_effects.csv',
            '.work/diagnostics/memory.csv',
        ],
        prepare=[
            'figure5/memory_scaling.csv',
            'figure5/memory_sample_flow.csv',
            'figure5/memory_parent_validation.csv',
            'private/figure5/memory_site_level_mapping.csv',
            'private/figure5/memory_model_frames.csv',
        ],
    ),
    'parameters': stage(
        None, 'figure4_trajectories', [],
        (), ['data.table'], None, None,
        deps=['plsc', 'tpv'],
        final=['private/figure4/primary_parameters.csv'],
    ),
}
groups = {name: [u for u, spec in STAGES.items() if spec["group"] == name]
          for name in dict.fromkeys(s["group"] for s in STAGES.values()) if name}
all_units = list(dict.fromkeys(u for g in groups.values() for u in g))
requested = x.stage if x.command in ["check", "_worker", "stage"] else x.command
if "," in requested:
    targets = requested.split(",")
    if (
        x.command not in ["check", "stage"]
        or len(targets) != len(set(targets))
        or any(u not in STAGES for u in targets)
    ):
        a.error("Stage list must contain distinct known calculation stages")
elif requested == "all":
    targets = all_units
elif requested in groups:
    targets = groups[requested]
elif requested in STAGES:
    targets = [requested]
else:
    a.error("Unknown stage")
units = []


def include(u):
    for d in STAGES[u]["deps"]:
        include(d)
    if u not in units:
        units.append(u)


for u in targets:
    include(u)
files = sorted(set(f for u in units for f in STAGES[u]["inputs"]))
contract = sorted(set(f for spec in STAGES.values() for f in spec["inputs"]))


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as s:
        for block in iter(lambda: s.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


RUNTIME_FILES = (
    "run.sh",
    "common.R",
    "figure1.R",
    "figure2.R",
    "figure3_associations.R",
    "figure3_lag.R",
    "figure4_trajectories.R",
    "figure4_plsc.m",
    "figure4_molecular.py",
    "figure5_progression.R",
    "figure5_mediation.R",
    "figure5_memory.R",
    "table_views.py",
    "matlab/rri_boot_check.m",
)


def runtime_sources():
    """Controllers and workers bind the same exact runtime file set."""
    found = {
        str(path.relative_to(code))
        for path in code.rglob("*")
        if path.is_file()
        and path.suffix in {".R", ".py", ".m", ".sh"}
        and "tests" not in path.relative_to(code).parts
    }
    expected = set(RUNTIME_FILES)
    if found != expected:
        missing = ", ".join(sorted(expected - found)) or "none"
        unexpected = ", ".join(sorted(found - expected)) or "none"
        raise ValueError(
            "Runtime source set changed; missing: " + missing + "; unexpected: " + unexpected
        )
    return {name: sha(code / name) for name in sorted(RUNTIME_FILES)}


def run(args, log=None):
    if log:
        with open(log, "w") as f:
            r = subprocess.run(args, stdout=f, stderr=subprocess.STDOUT)
    else:
        r = subprocess.run(args)
    if r.returncode:
        raise RuntimeError("Command failed (see local log): " + shlex.join(args[:2]))


def capture(args):
    return subprocess.check_output(args, text=True, stderr=subprocess.PIPE).strip()


def r_packages(names):
    expr = (
        "p<-c(" + ",".join(repr(p) for p in names) + ");"
        'cat(as.character(getRversion()),"\\n");'
        'for(n in p){f<-tryCatch(find.package(n),error=function(e)"");'
        'cat(n,if(nzchar(f))as.character(packageVersion(n))else"ABSENT",f,sep="\\t");'
        'cat("\\n")}'
    )
    return capture([rscript, "--vanilla", "-e", expr])


def package_artifacts(root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*"))
            if p.is_file() and (p.suffix in {".rdb", ".rdx", ".so"}
                               or p.name in {"DESCRIPTION", "NAMESPACE"})}


def verify_installation():
    """Verify the startup software/source manifest again before releasing outputs."""
    if runtime_sources() != identity["source"]:
        raise ValueError("Scientific source changed during execution")
    for key, digest_key in (("rscript", "rscript_sha256"),
                            ("matlab", "matlab_executable_sha256")):
        if runtime.get(key) and sha(runtime[key]) != runtime[digest_key]:
            raise ValueError("Numerical executable changed during execution: " + key)
    if runtime.get("R_packages"):
        names = [row.split("\t")[0] for row in runtime["R_packages"].splitlines()[1:]]
        if r_packages(names) != runtime["R_packages"]:
            raise ValueError("R package version/path changed during execution")
        for row in runtime["R_packages"].splitlines()[1:]:
            fields = row.split("\t")
            if len(fields) >= 3 and fields[2]:
                if package_artifacts(P(fields[2])) != runtime["R_package_artifacts"][fields[0]]:
                    raise ValueError("R installation changed during execution: " + fields[0])
    if runtime.get("SEM_artifacts") is not None:
        root = P(runtime["SEM_LIBRARY"]) / "lavaan"
        if package_artifacts(root) != runtime["SEM_artifacts"]:
            raise ValueError("SEM installation changed during execution")
    if toolbox:
        root = P(toolbox).resolve(strict=True)
        observed = {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*.m"))}
        if observed != runtime["toolbox"]:
            raise ValueError("PLS toolbox changed during execution")
    import importlib.metadata as md
    for package, expected in runtime["numeric_python"].items():
        try:
            actual = md.version(package)
        except md.PackageNotFoundError:
            actual = "ABSENT"
        if actual != expected:
            raise ValueError("Python package changed during execution: " + package)


def atomic_json(path, value):
    t = path.with_suffix(path.suffix + ".tmp")
    t.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")
    t.replace(path)


# Controllers bind inputs and software; workers bind the prepared objects they consume.
if x.command == "_worker":
    if x.stage not in STAGES or not x.task or x.task < 1:
        a.error("Internal worker needs stage/task")
    identity = json.loads((work / "identity.json").read_text())
    source = runtime_sources()
    if source != identity["source"]:
        raise ValueError("Worker scientific source identity changed")
    runtime = identity["runtime"]
    for key, value in {
        "python": sys.version,
        "python_executable": sys.executable,
        "rscript": shutil.which(rscript),
        "R_LIBS_USER": os.environ.get("R_LIBS_USER", ""),
        "R_LIBS_SITE": os.environ.get("R_LIBS_SITE", ""),
        "SEM_LIBRARY": os.environ.get("FIGURE5_SEM_LIBRARY", ""),
    }.items():
        if runtime.get(key) != value:
            raise ValueError("Worker runtime setting changed: " + key)
    if x.stage not in ["plsc", "molecular"]:
        if sha(runtime["rscript"]) != runtime["rscript_sha256"]:
            raise ValueError("Worker R executable changed")
        selected = set(STAGES[x.stage]["packages"]) - {"lavaan"}
        if selected & {"lme4", "lmerTest", "mgcv"}:
            selected.update(["Matrix", "MASS", "reformulas"])
        observed = r_packages(sorted(selected)).splitlines()
        expected = runtime["R_packages"].splitlines()
        expected_rows = {z.split("\t")[0]: z for z in expected[1:]}
        if observed[0].strip() != expected[0].strip() or [
            z.rstrip("\t") for z in observed[1:]
        ] != [expected_rows[n].rstrip("\t") for n in sorted(selected)]:
            raise ValueError("Worker stage package version/path changed")
        # The SEM task verifies its explicitly resolved namespace and package bytes in-process.
    if x.stage == "plsc":
        if (
            shutil.which(matlab) != runtime["matlab"]
            or sha(runtime["matlab"]) != runtime["matlab_executable_sha256"]
        ):
            raise ValueError("Worker MATLAB executable changed")
        for name, digest in runtime["toolbox"].items():
            if not (P(toolbox) / name).is_file() or sha(P(toolbox) / name) != digest:
                raise ValueError("Worker toolbox changed: " + name)
    if x.stage == "molecular":
        import importlib.metadata as md

        if {n: md.version(n) for n in ["numpy", "pandas", "scipy"]} != runtime["numeric_python"]:
            raise ValueError("Worker molecular runtime changed")
    signature = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    if (out / "run_signature.txt").read_text().strip() != signature:
        raise ValueError("Worker controller identity receipt changed")
    os.environ["MODEL_SIGNATURE"] = signature
else:
    inputs = {}
    for name in contract:
        p = data / name
        if p.exists() or p.is_symlink():
            if not p.is_file() or not os.access(p, os.R_OK):
                raise ValueError("Unreadable/invalid input " + name)
            with p.open(newline="") as h:
                reader = csv.reader(h)
                header = next(reader, None)
                if not header or len(header) != len(set(header)):
                    raise ValueError("Invalid header " + name)
                rows = 0
                for row in reader:
                    if len(row) != len(header):
                        raise ValueError("Ragged input " + name)
                    rows += 1
                if not rows:
                    raise ValueError("Empty supported input " + name)
            inputs[name] = {"sha256": sha(p), "rows": rows, "columns": header}
        else:
            inputs[name] = {"absent": True}
    for name in files:
        if inputs[name].get("absent"):
            raise ValueError("Missing required input " + name)
    if (data / "input_manifest.csv").is_file():
        with (data / "input_manifest.csv").open(newline="") as h:
            manifest = list(csv.DictReader(h))
        if len(manifest) != len({r["relative_file"] for r in manifest}):
            raise ValueError("Duplicate input manifest file")
        for r in manifest:
            name = r["relative_file"]
            v = inputs.get(name)
            if (
                v is None
                or v.get("absent")
                or v["sha256"] != r["sha256"]
                or v["rows"] != int(r["rows"])
                or len(v["columns"]) != int(r["columns"])
            ):
                raise ValueError("Input manifest mismatch " + name)
    # Bind additions as well as deletions. Data documentation changes also invalidate resume.
    extra = {p.name: sha(p) for p in sorted(data.glob("*.csv")) if p.name not in contract}
    source = runtime_sources()
    packages = [
        "data.table",
        "lme4",
        "lmerTest",
        "mgcv",
        "MASS",
        "Matrix",
        "digest",
        "jsonlite",
        "lavaan",
        "reformulas",
    ]
    # Missing optional engines are recorded, not required for unrelated figures.
    runtime = {
        "python": sys.version,
        "python_executable": sys.executable,
        "rscript": shutil.which(rscript),
        "matlab": shutil.which(matlab),
        "R_LIBS_USER": os.environ.get("R_LIBS_USER", ""),
        "SEM_LIBRARY": os.environ.get("FIGURE5_SEM_LIBRARY", ""),
        "R_LIBS_SITE": os.environ.get("R_LIBS_SITE", ""),
    }
    if runtime["rscript"]:
        runtime["R_packages"] = r_packages(packages)
        runtime["rscript_sha256"] = sha(runtime["rscript"])
    if toolbox:
        tb = P(toolbox).resolve(strict=True)
        runtime["toolbox"] = {str(p.relative_to(tb)): sha(p) for p in sorted(tb.rglob("*.m"))}
    if runtime["matlab"]:
        runtime["matlab_executable_sha256"] = sha(runtime["matlab"])
    if any(u not in ["plsc", "molecular"] for u in units) and not runtime["rscript"]:
        raise ValueError("Set RSCRIPT_BIN to existing R4.4.1")

    if runtime["rscript"]:
        probe = subprocess.run(
            [rscript, "--vanilla", str(code / "figure5_mediation.R"), "--action", "runtime"],
            text=True,
            capture_output=True,
        )
        runtime["SEM_resolved"] = (
            json.loads(probe.stdout)
            if probe.returncode == 0
            else {"available": False, "reason": probe.stderr.strip()}
        )
    else:
        runtime["SEM_resolved"] = {"available": False, "reason": "Rscript unavailable"}
    if "mediation" in units and not runtime["SEM_resolved"].get("available"):
        raise ValueError("SEM preflight failed: " + runtime["SEM_resolved"]["reason"])

    versions = {
        z[0]: z[1]
        for line in runtime.get("R_packages", "").splitlines()[1:]
        for z in [line.split("\t")]
        if len(z) >= 2
    }
    if runtime["SEM_resolved"].get("available"):
        versions["lavaan"] = runtime["SEM_resolved"]["version"]
    for package in sorted({p for u in units for p in STAGES[u]["packages"]}):
        if versions.get(package, "ABSENT") == "ABSENT":
            raise ValueError("Missing stage-required R package: " + package)
    if (
        runtime.get("R_packages", "").splitlines()
        and runtime["R_packages"].splitlines()[0].strip() != "4.4.1"
    ):
        raise ValueError("This reviewed runtime requires R 4.4.1")
    if "plsc" in units and (not runtime["matlab"] or not toolbox):
        raise ValueError(
            "Set MATLAB_BIN and PLS_TOOLBOX_ROOT to existing licensed MATLAB R2019b/original toolbox"
        )
    import importlib.metadata as md

    runtime["numeric_python"] = {}
    for package in ["numpy", "pandas", "scipy"]:
        try:
            runtime["numeric_python"][package] = md.version(package)
        except md.PackageNotFoundError:
            runtime["numeric_python"][package] = "ABSENT"
    if "molecular" in units and "ABSENT" in runtime["numeric_python"].values():
        raise ValueError("Existing numpy, pandas and scipy required")
    # Hash installed numerical code once at startup and again before releasing outputs.
    runtime["R_package_artifacts"] = {}
    if runtime.get("R_packages"):
        for line in runtime["R_packages"].splitlines()[1:]:
            fields = line.split("\t")
            if len(fields) >= 3 and fields[2]:
                root = P(fields[2])
                runtime["R_package_artifacts"][fields[0]] = package_artifacts(root)
    sem = os.environ.get("FIGURE5_SEM_LIBRARY", "")
    if sem:
        runtime["SEM_artifacts"] = package_artifacts(P(sem) / "lavaan")
    identity = {
        "schema": "main_csv_replication_v3",
        "inputs": inputs,
        "extra_inputs": extra,
        "source": source,
        "runtime": runtime,
        "profile": "analysis",
        "scope": "main",
        "scientific_counts": {
            "ICC_bootstraps": 5000,
            "PLSC_permutations": 5000,
            "PLSC_bootstraps": 5000,
            "GAM_draws": 10000,
            "Moran_draws": 10000,
        },
    }
    signature = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    identity_file = work / "identity.json"
    if identity_file.exists():
        if json.loads(identity_file.read_text()) != identity:
            raise ValueError(
                "Input/code/runtime/settings identity changed; use a fresh output root"
            )
    else:
        if any(work.glob("*/prepared.rds")):
            raise ValueError("Unsigned existing work is not reusable")
        atomic_json(identity_file, identity)
    os.environ["MODEL_SIGNATURE"] = signature
    signature_path = out / "run_signature.txt"
    if signature_path.exists():
        if signature_path.read_text().strip() != signature:
            raise ValueError("Root signature differs")
    else:
        signature_path.write_text(signature + "\n")
# Same-run dependencies cannot inherit a signature override from the shell.
os.environ.pop("MODEL_TRAJECTORY_PRODUCER_SIGNATURE", None)
os.environ.pop("MODEL_TPV_PRODUCER_SIGNATURE", None)
if x.command == "check":
    print(f"CHECK PASS; scope={requested}; required inputs={len(files)}; no study model executed")
    print("Stages: " + ", ".join(units))
    sys.exit(0)


def command(unit, action, task=None):
    spec = STAGES[unit]
    if spec["engine"] == "MATLAB":
        q = lambda s: "'" + str(s).replace("'", "''") + "'"
        call = (
            "addpath("
            + q(code)
            + ");" + spec["script"] + "("
            + ",".join([q(action), q(out), q(toolbox)] + ([str(task)] if task else []))
            + ");"
        )
        return [matlab, "-singleCompThread", "-batch", call]
    if spec["engine"] == "Python":
        c = [
            sys.executable,
            str(code / (spec["script"] + ".py")),
            "--mode",
            action,
            "--data",
            str(data),
            "--out",
            str(out),
        ]
        return c + (["--task-id", str(task)] if task else [])
    c = [rscript, "--vanilla", str(code / (spec["script"] + ".R")), "--action", action]
    c += list(spec["flags"])
    if task:
        c += ["--task-id" if unit in ["trajectories", "tpv"] else "--task", str(task)]
    return c


def completed_receipt(unit, task=None):
    return work / "receipts" / (unit + ("_" + str(task) if task else "") + ".json")


class MissingCompletionFile(FileNotFoundError):
    """A declared completion file is not visible yet."""

    def __init__(self, message, paths):
        self.paths = tuple(map(str, paths))
        super().__init__(2, message, self.paths[0])


# Output ownership is explicit, including files rewritten by interrupted collectors.
def csv_rows(path, retain=True):
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or any(not c for c in header) or len(header) != len(set(header)):
            raise ValueError("Missing/invalid output schema " + str(path))
        rows, count = [], 0
        for row in reader:
            if len(row) != len(header):
                raise ValueError("Ragged stage output " + str(path))
            count += 1
            if retain:
                rows.append(dict(zip(header, row)))
    return header, rows if retain else count


def task_ids(unit):
    if unit == "molecular":
        ids = json.loads((work / unit / "prepared.json").read_text())["selected_task_ids"]
        if ids != list(STAGES[unit]["tasks"]):
            raise ValueError("Incomplete main molecular task plan")
        return ids
    _, rows = csv_rows(work / unit / "tasks.csv")
    ids = [int(r["task_id"]) for r in rows]
    if ids != list(STAGES[unit]["tasks"]) or int(
        (work / unit / "n_tasks.txt").read_text()
    ) != len(ids):
        raise ValueError("Incomplete/invalid task plan " + unit)
    return ids


def owned_outputs(unit, action, task=None):
    absent = []
    if action == "prepare":
        names = list(STAGES[unit]["prepare"])
        base = ".work/" + unit + "/"
        if unit == "plsc":
            names += [
                base + s
                for s in [
                    "prepared_manifest.mat",
                    "prepared_manifest.json",
                    "tasks.csv",
                    "n_tasks.txt",
                ]
            ] + [
                base + "inputs/" + n + "_lme_" + p + ".mat"
                for n in ["LC", "SNVTA"]
                for p in ["intercept", "slope"]
            ]
        elif unit == "molecular":
            names += [
                base + s
                for s in [
                    "prepared.json",
                    "preparation_receipt.json",
                    "n_tasks.txt",
                    "molecular_observed.csv",
                    "molecular_unestimated.csv",
                    "molecular_plot_coordinates.csv",
                    "molecular_tasks.csv",
                    "molecular_null_chunks.csv",
                ]
            ]
            prepared = json.loads((work / unit / "prepared.json").read_text())
            names += [
                str(P(name).resolve().relative_to(out)) for name in prepared["prepared_identity"]
            ]
        else:
            names += [base + s for s in ["prepared.rds", "tasks.csv", "n_tasks.txt"]]
        if unit in ["concurrent", "prospective", "lag"]:
            names += [base + f"outcome_{j:02d}.rds" for j in range(1, 39)]
        task_ids(unit)
    elif action == "task":
        if task not in task_ids(unit):
            raise ValueError("Unsupported task identity " + unit)
        base = ".work/" + unit + "/"
        if unit == "plsc":
            names = [base + f"task_{task:04d}.json"] + [
                base + f"tasks/{task:04d}/" + s
                for s in [
                    "result.mat",
                    "global.csv",
                    "regional.csv",
                    "scores.csv",
                    "behavior.csv",
                    "runtime_receipt.json",
                    "score_transformation.json",
                    "transforms.csv",
                ]
            ]
        elif unit == "molecular":
            name = base + f"task_{task:04d}.json"
            names = [name]
            if not os.path.lexists(out / name):
                raise MissingCompletionFile("Missing molecular task record", [out / name])
            d = json.loads((out / name).read_text())
            if d.get("task_id") != task or d.get("status") != "complete":
                raise ValueError("Molecular task is incomplete")
            if d.get("model_status") == "estimated":
                names += [
                    base + f"tasks/task_{task:03d}.json",
                    base + f"tasks/task_{task:03d}.npz",
                ]
            elif d.get("model_status") != "not_estimable" or not d.get("reason"):
                raise ValueError("Molecular unavailable task lacks a scientific reason")
        else:
            names = [base + f"task_{task:04d}.rds"]
    else:
        names = list(STAGES[unit]["final"])
        if unit in ["concurrent", "prospective", "lag"]:
            # Bind both present and scientifically absent model-frame exports.
            _, diagnostics = csv_rows(out / (".work/diagnostics/figure3_" + unit + ".csv"))
            has_estimate = any(r.get("status") == "estimated" for r in diagnostics)
            for name in [
                "private/figure3/" + unit + "_model_frames.csv",
                "private/figure3/" + unit + "_factor_references.csv",
                "figure3/" + unit + "_scale_references.csv",
            ]:
                (names if has_estimate or (out / name).is_file() else absent).append(name)
    return sorted(set(names)), sorted(absent)


def validate_coverage(unit):
    def keys(name, columns, expected):
        header, rows = csv_rows(out / name)
        if not set(columns) <= set(header):
            raise ValueError("Missing coverage columns " + name)
        values = [tuple(r[c] for c in columns) for r in rows]
        if len(values) != len(expected) or set(values) != set(expected):
            raise ValueError("Incomplete/duplicate scientific output coverage " + name)

    def groups(name, columns, expected):
        header, rows = csv_rows(out / name)
        if not set(columns) <= set(header) or {tuple(r[c] for c in columns) for r in rows} != set(
            expected
        ):
            raise ValueError("Incomplete model status coverage " + name)

    if unit == "technical":
        groups(
            ".work/diagnostics/technical.csv",
            ["nucleus", "field_strength", "hemisphere"],
            [(n, f, h) for n in ["LC", "SNVTA"] for f in ["1.5T", "3T"] for h in ["L", "R"]],
        )
    if unit == "figure2":
        groups(".work/diagnostics/figure2.csv", ["nucleus"], [("LC",), ("SNVTA",)])
    if unit == "icc":
        keys(
            "figure1/icc.csv",
            ["nucleus", "field_strength", "hemisphere"],
            [(n, f, h) for n in ["LC", "SNVTA"] for f in ["1.5T", "3T"] for h in ["L", "R"]],
        )
    if unit in ["concurrent", "prospective", "lag"]:
        _, plan = csv_rows(work / unit / "tasks.csv")
        expected = [(r["outcome"], r["nucleus"]) for r in plan]
        keys(
            "figure3/" + unit + ("_models.csv" if unit == "lag" else "_associations.csv"),
            ["outcome", "nucleus"],
            expected,
        )
        if unit == "lag":
            keys("figure3/lag_interval_summary.csv", ["outcome", "nucleus"], expected)
    if unit == "plsc":
        keys("figure4/plsc_global_results.csv", ["task_id"], [(str(i),) for i in STAGES[unit]["tasks"]])
        keys(
            "figure4/plsc_regional_results.csv",
            ["task_id", "region_index"],
            [(str(i), str(j)) for i in STAGES[unit]["tasks"] for j in range(1, 146)],
        )
    if unit == "molecular":
        _, pet = csv_rows(data / "pet_regions.csv")
        ids = sorted({r["pet_id"] for r in pet})
        if len(ids) != 22:
            raise ValueError("Expected 22 frozen PET map identities")
        keys(
            "figure4/molecular_results.csv",
            ["branch", "pet_id"],
            [(n, p) for n in ["LC", "SNVTA"] for p in ids],
        )
    if unit == "mediation":
        keys(
            "figure5/mediation_effects.csv",
            ["task_id", "model_variant", "effect"],
            [
                (str(i), v, e)
                for i in STAGES[unit]["tasks"]
                for v in ["same_frame_base", "same_frame_tpv"]
                for e in ["indirect", "direct", "total"]
            ],
        )
    if unit == "amyloid":
        keys(
            "figure5/memory_focal_effects.csv",
            ["task_id", "effect_role"],
            [
                (str(i), e)
                for i in STAGES[unit]["tasks"]
                for e in ["local_pathology_time", "plsc_pathology_time", "tpv_pathology_time"]
            ],
        )


def output_state(unit, action, task=None, expected=None):
    names, absent = owned_outputs(unit, action, task)
    outputs = {}
    schemas = {}
    if not names:
        raise ValueError("Stage declares no owned output " + unit)
    if expected is not None and (
        set(expected["outputs"]) != set(names)
        or set(expected["schemas"]) != {name for name in names if P(name).suffix == ".csv"}
        or expected.get("absent_outputs") != absent
    ):
        raise ValueError("Receipt owned output changed/missing " + unit)
    missing = []
    for name in names:
        p = out / name
        if not os.path.lexists(p):
            missing.append(p)
            continue
        if not p.is_file() or not os.access(p, os.R_OK) or p.stat().st_size == 0:
            raise ValueError("Invalid/unreadable owned output " + name)
        if p.suffix == ".csv":
            header, count = csv_rows(p, retain=False)
            schemas[name] = {"columns": header, "rows": count}
        outputs[name] = sha(p)
        if expected is not None and (
            expected["outputs"][name] != outputs[name]
            or (name in schemas and expected["schemas"][name] != schemas[name])
        ):
            raise ValueError("Receipt owned output changed/missing " + unit + ": " + name)
    if missing:
        raise MissingCompletionFile("Missing owned output", missing)
    if action not in ["prepare", "task"]:
        validate_coverage(unit)
    return outputs, schemas, absent


def verify_receipt(unit, action, task=None):
    receipt = (
        completed_receipt(unit + "_prepare")
        if action == "prepare"
        else completed_receipt(unit, task)
    )
    if not os.path.lexists(receipt):
        raise MissingCompletionFile("Missing completion receipt", [receipt])
    if not receipt.is_file():
        raise ValueError("Invalid completion receipt " + receipt.name)
    d = json.loads(receipt.read_text())
    if (
        d.get("signature") != signature
        or d.get("unit") != unit
        or d.get("task") != task
        or d.get("action") != action
        or d.get("exit_code") != 0
    ):
        raise ValueError("Worker receipt identity/completion changed")
    if (
        not isinstance(d.get("outputs"), dict)
        or not d["outputs"]
        or any(not isinstance(v, str) or len(v) != 64 or any(c not in "0123456789abcdef" for c in v)
               for v in d["outputs"].values())
        or not isinstance(d.get("schemas"), dict)
        or not isinstance(d.get("absent_outputs"), list)
    ):
        raise ValueError("Malformed receipt output metadata " + unit)
    for schema in d["schemas"].values():
        if not isinstance(schema, dict) or set(schema) != {"columns", "rows"}:
            raise ValueError("Malformed receipt schema " + unit)
        columns, rows = schema["columns"], schema["rows"]
        if (
            not isinstance(columns, list) or not columns
            or any(not isinstance(c, str) or not c for c in columns)
            or len(columns) != len(set(columns))
            or type(rows) is not int or rows < 0
        ):
            raise ValueError("Malformed receipt schema " + unit)
    output_state(unit, action, task, expected=d)


def verify_scheduler_completion(unit, ids, timeout=None, delay=2):
    # Recheck only absent files; every arrival still requires full verification.
    if timeout is None:
        timeout = completion_wait
    deadline = time.monotonic() + timeout
    pending = dict.fromkeys(ids)
    while pending:
        for task, missing in list(pending.items()):
            if missing is not None and not any(os.path.lexists(path) for path in missing):
                continue
            checking_started = time.monotonic()
            try:
                verify_receipt(unit, "task", task)
            except MissingCompletionFile as error:
                pending[task] = error.paths
            else:
                del pending[task]
            finally:
                # Integrity checks have no visibility-wait deadline.
                deadline += time.monotonic() - checking_started
        if not pending:
            return
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise RuntimeError(
                "Completion visibility timeout; no automatic resubmission: "
                + "; ".join(f"{unit} task {task}: {', '.join(paths)}" for task, paths in pending.items())
            )
        time.sleep(min(delay, remaining))


def invoke(unit, action, task=None, tasks_verified=False):
    receipt = (
        completed_receipt(unit + "_prepare")
        if action == "prepare"
        else completed_receipt(unit, task)
    )
    if action == "finalize" and not tasks_verified:
        for i in task_ids(unit):
            verify_receipt(unit, "task", i)
    if receipt.exists():
        verify_receipt(unit, action, task)
        return
    log = work / "logs" / (unit + "_" + action + ("_" + str(task) if task else "") + ".log")
    started = time.monotonic()
    run(command(unit, action, task), log)
    outputs, schemas, absent = output_state(unit, action, task)
    atomic_json(
        receipt,
        {
            "signature": signature,
            "unit": unit,
            "action": action,
            "task": task,
            "seconds": time.monotonic() - started,
            "exit_code": 0,
            "outputs": outputs,
            "schemas": schemas,
            "absent_outputs": absent,
        },
    )


def submit_tasks(unit, ids):
    if not shutil.which("qsub"):
        raise RuntimeError("SGE requested but qsub unavailable")
    pending = [i for i in ids if not completed_receipt(unit, i).exists()]
    # Submit only missing tasks. -sync waits for scheduler return; receipts prove completion.
    if pending:
        worker = work / (unit + "_worker.sh")
        exports = [
            "RSCRIPT_BIN",
            "R_LIBS_USER",
            "R_LIBS_SITE",
            "FIGURE5_SEM_LIBRARY",
            "MATLAB_BIN",
            "PLS_TOOLBOX_ROOT",
            "PYTHON_BIN",
        ]
        body = "#!/bin/bash\nset -euo pipefail\n" + "".join(
            "export " + n + "=" + shlex.quote(os.environ[n]) + "\n"
            for n in exports
            if n in os.environ
        )
        body += "ids=(" + " ".join(map(str, pending)) + ")\n"
        body += (
            "exec bash "
            + shlex.quote(str(code / "run.sh"))
            + " _worker --data "
            + shlex.quote(str(data))
            + " --out "
            + shlex.quote(str(run_root))
            + " --stage "
            + unit
            + ' --task "${ids[$SGE_TASK_ID-1]}"\n'
        )
        worker.write_text(body)
        worker.chmod(0o700)
        memory_default, hours_default = STAGES[unit]["resources"]
        memory = os.environ.get("SGE_MEMORY", memory_default)
        hours = os.environ.get("SGE_WALLTIME", hours_default)
        if not memory.endswith("G") or not 0 < float(memory[:-1]) <= 24:
            raise ValueError("SGE_MEMORY must be <=24G")
        h, m, s = map(int, hours.split(":"))
        if h * 3600 + m * 60 + s > 43200:
            raise ValueError("SGE_WALLTIME must be <=12h")
        task_limit = min(x.jobs, 2) if STAGES[unit]["engine"] == "MATLAB" else x.jobs
        cmd = ["qsub", "-sync", "y", "-cwd", "-V",
               "-t", "1-" + str(len(pending)), "-tc", str(task_limit),
               "-pe", os.environ.get("SGE_PE", "smp"), "1",
               "-l", "h_vmem=" + memory + ",h_rt=" + hours,
               "-o", str(work / "logs"), "-e", str(work / "logs")]

        if os.environ.get("SGE_QUEUE"):
            cmd += ["-q", os.environ["SGE_QUEUE"]]
        run(cmd + [str(worker)], work / "logs" / (unit + "_qsub.log"))
    verify_scheduler_completion(unit, ids)


if x.command == "_worker":
    verify_receipt(x.stage, "prepare")
    invoke(x.stage, "task", x.task)
    sys.exit(0)
# Serial stage dependencies; concurrency only distributes whole existing tasks.
for unit in units:
    if unit == "parameters":
        invoke(unit, "primary-parameters")
        continue
    if completed_receipt(unit).exists():
        invoke(unit, "prepare")
        invoke(unit, "finalize")
        continue
    invoke(unit, "prepare")
    ids = task_ids(unit)
    if x.mode == "local":
        for task in ids:
            invoke(unit, "task", task)
    else:
        submit_tasks(unit, ids)
    invoke(unit, "finalize", tasks_verified=True)
    print(f"{unit}: {len(ids)} model tasks completed.")
verify_installation()
run(
    [rscript, "--vanilla", str(code / "common.R"), "--action", "diagnostics"],
    work / "logs" / "diagnostics.log",
)
sys.path.insert(0, str(code))
from table_views import export_views

views = export_views(
    run_root / "results",
    source_root=out,
    stages=units,
    receipt=work / "receipts" / ("table_views_" + requested.replace(",", "_") + ".json"),
)
print(f"Results written to {run_root / 'results'} ({len(views['files'])} CSVs).")
PY
