#!/usr/bin/env python3
"""One analysis graph, executed locally or through native sbatch. Plans by default."""
import argparse
import csv
import contextlib
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import uuid

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MODES = ("main", "no_ces", "no_overlap", "age")
# Readers depend explicitly on every producer. Independent summaries may overlap.
PARENTS = {
    "prepare": [], "model_fit": ["prepare"], "models": ["model_fit"],
    "simulation_prepare": ["prepare"], "simulation_fit": ["simulation_prepare"],
    "simulations": ["simulation_fit"], "upper_prepare": ["prepare"],
    "upper_fit": ["upper_prepare"], "upper": ["upper_fit"],
    "korean": ["prepare"], "public_example": ["prepare"],
    "main": ["models"], "cohort": ["models"], "ddid": ["models"],
    "secondary": ["models"], "posterior": ["models"], "prevalence": ["models"],
    "forecast": ["models"], "dd": ["main"], "no_overlap": ["dd"],
    "paternal": ["cohort"], "figure2": ["main", "cohort"],
    "figure3": ["main", "dd", "forecast"], "figure4": ["dd"],
    "supplementary": ["secondary", "posterior", "ddid", "prevalence", "no_overlap", "paternal", "korean"],
    "tables": ["main", "cohort", "dd", "forecast", "simulations"],
    "checklist": ["tables", "posterior", "secondary", "paternal", "upper", "ddid"],
    "finalize": ["figure2", "figure3", "figure4", "supplementary", "checklist", "public_example"],
}
ARRAYS = {"model_fit": 245, "simulation_fit": 500, "upper_fit": 30}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{uuid.uuid4().hex}.tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)


def code_hashes():
    paths = []
    for folder in ("R", "analysis", "inst/examples"):
        paths += [p for p in (REPO / folder).rglob("*")
                  if p.is_file() and p.suffix in (".R", ".py", ".sh", ".json", ".tsv")
                  and not any(x in p.parts for x in ("__pycache__", "inputs", "outputs", "logs"))]
    return {str(p.relative_to(REPO)): sha(p) for p in sorted(paths)}


def load_config(path):
    c = json.loads(Path(path).read_text())
    if Path(path).resolve() == (HERE / "manuscript.json").resolve():
        for key, name in {"root": "BURDENMLEDN_ANALYSIS_ROOT", "run_id": "RUN_DATE",
                          "bootstrap_draws": "BURDENMLEDN_PREVIOUS_BOOTSTRAPS", "rscript": "R_BIN"}.items():
            if os.environ.get(name): c[key] = os.environ[name]
    if c["profile"] != "manuscript-r2" or not re.fullmatch(r"[A-Za-z0-9_-]+", c["run_id"]):
        raise ValueError("Invalid profile/run identifier")
    # This profile represents the approved study, not an arbitrary experiment.
    fixed = dict(bootstraps=100, simulation_replicates=500, comparison_replicates=100,
                 simulation_batch_size=25, asd_rr_bound=100, dd_rr_bound=250)
    for key, value in fixed.items():
        if c[key] != value:
            raise ValueError(f"manuscript-r2 requires {key}={value}")
    if not c["root"] or not Path(c["root"]).is_absolute():
        raise ValueError("Set root to an absolute analysis directory")
    root = Path(c["root"]).resolve()
    if root == REPO or REPO in root.parents:
        raise ValueError("Keep run outputs outside the source checkout")
    if c["bootstrap_draws"] and not Path(c["bootstrap_draws"]).is_absolute():
        raise ValueError("bootstrap_draws must be an absolute path or empty")
    for key in ("model_seed", "simulation_seed", "korean_seed"):
        if type(c[key]) is not int or not 0 <= c[key] < 2**31 - 501:
            raise ValueError("Invalid integer seed: " + key)
    if type(c["forecast_cores"]) is not int or c["forecast_cores"] < 1:
        raise ValueError("forecast_cores must be positive")
    return c


def environment(c):
    env = os.environ.copy()
    root = Path(c["root"])
    env.update(BURDENMLEDN_ANALYSIS_ROOT=str(root), BURDENMLEDN_REPO_ROOT=str(REPO),
               RUN_DATE=c["run_id"], MODEL_SEED=str(c["model_seed"]),
               SIMULATION_SEED=str(c["simulation_seed"]),
               BURDENMLEDN_PREVIOUS_BOOTSTRAPS=c["bootstrap_draws"],
               R_LIBS=str(root / "R-library") + (os.pathsep + env["R_LIBS"] if env.get("R_LIBS") else ""),
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    return env


def commands(stage, c, task=None):
    root = Path(c["root"])
    scripts = REPO / "analysis/scripts"
    r = lambda name, *args: [c["rscript"], "--vanilla", str(scripts / (name + ".R")), *map(str, args)]
    w = lambda name, *args: [c["rscript"], "--vanilla", str(HERE / (name + ".R")), *map(str, args)]
    manifest = lambda mode: str(root / f"outputs/data/model_manifest_{mode}_{c['run_id']}.rds")
    model = lambda mode: r("run_models", "--mode", mode, "--run-date", c["run_id"],
                           "--seed", c["model_seed"], "--bootstraps", c["bootstraps"],
                           "--run-autism", "true", "--run-ddd", "true",
                           "--autism-max-effect-size", c["asd_rr_bound"],
                           "--ddd-max-effect-size", c["dd_rr_bound"])
    if stage == "prepare":
        return [w("preflight")] + [model(m) for m in MODES]
    if stage == "model_fit":
        paths = []
        for mode in MODES:
            paths += (root / f"jobs/{c['run_id']}/{mode}/tasks.txt").read_text().splitlines()
        if len(paths) != ARRAYS[stage]:
            raise ValueError("Model inventory does not match the manuscript profile")
        return [r("run_model_job", paths[task - 1])]
    if stage == "models":
        return [model(m) for m in MODES]
    if stage.startswith("simulation_"):
        return [w("simulation_batches", "prepare" if stage.endswith("prepare") else "worker", *([] if task is None else [task]))]
    if stage == "simulations":
        return [w("simulation_batches", "collect"),
                r("analyze_simulations", "--run-date", c["run_id"], "--seed", c["simulation_seed"],
                  "--optimizer", "mixsqp", "--compare-optimizer", "EM", "--comparison-replicates", 100)]
    if stage.startswith("upper_") or stage == "upper":
        action = {"upper_prepare": "prepare", "upper_fit": "worker", "upper": "collect"}[stage]
        return [w("upper_bound_jobs", action, *([] if task is None else [task]))]
    if stage == "korean":
        return [r("run_korean_wgs", "--run-date", c["run_id"], "--seed", c["korean_seed"], "--n-boot", 100)]
    if stage == "public_example":
        return [[c["rscript"], "--vanilla", str(REPO / "inst/examples/kaplanis/run_kaplanis_public_example.R")]]
    summary = {"main": "main_autism", "cohort": "cohort_sex", "ddid": "ddid",
               "secondary": "secondary_autism", "posterior": "posterior_diagnostics",
               "prevalence": "autism_dd_prevalence", "dd": "autism_dd"}
    if stage in summary:
        return [r("summarize_" + summary[stage], "--model-manifest", manifest("main"))]
    if stage == "no_overlap":
        return [r("summarize_no_overlap", "--model-manifest", manifest("no_overlap"))]
    if stage == "forecast":
        return [r("forecasting_script_revision", "--model-manifest", manifest("main"),
                  "--seed", c["model_seed"], "--cores", c["forecast_cores"])]
    if stage == "paternal":
        return [r("describe_us_father_age", "--cohort-summary-file", root / "outputs/derived/cohort_ddid/cohort_sex_summary.rds"),
                r("compare_spark_us_father_age")]
    if stage in ("figure2", "figure3", "figure4"):
        return [r({"figure2": "make_figure2_cohort_sex", "figure3": "make_forecasting_figures",
                   "figure4": "make_figure4_autism_dd"}[stage])]
    if stage == "supplementary":
        return [r("make_supplementary_figures", "--no-ces-model-manifest", manifest("no_ces"))]
    if stage in ("tables", "checklist"):
        return [r("make_supplementary_tables" if stage == "tables" else "make_manuscript_estimate_table",
                  "--model-manifest", manifest("main"))]
    if stage == "finalize":
        return [w("assemble_submission"), w("validate_submission")]
    raise ValueError(stage)


def products(stage, c, task=None):
    root = Path(c["root"])
    if stage == "model_fit":
        path = Path(commands(stage, c, task)[0][-1])
        return [path.with_name(path.stem + "_result.rds")]
    patterns = {
        "prepare": [f"jobs/{c['run_id']}/*/*.rds", f"jobs/{c['run_id']}/*/tasks.txt", "provenance/input_manifest.tsv"],
        "models": [f"outputs/data/models*_{c['run_id']}.Rdata", f"outputs/data/model_manifest_*_{c['run_id']}.rds"],
        "simulation_prepare": ["jobs/simulations/inputs.rds", "jobs/simulations/inventory.tsv", "outputs/data/effect_size_samples_*"],
        "simulation_fit": [f"jobs/simulations/batch_{task:03d}.rds"] if task else [],
        "simulations": ["outputs/tables/diagnostics/simulations/*", "outputs/tables/supplementary/SupplementaryTable1_*", "outputs/figures/main/Figure1.pdf", "outputs/figures/supplementary/SupplementaryFigure1_*", "outputs/figures/supplementary/SupplementaryFigure16_*"],
        "upper_prepare": ["outputs/derived/upper_rr_no_CES/workspace_private.rds", "outputs/derived/upper_rr_no_CES/inventory.tsv", "outputs/derived/upper_rr_no_CES/excluded_genes.tsv"],
        "upper_fit": [f"outputs/derived/upper_rr_no_CES/job_{task:02d}.rds"] if task else [],
        "upper": ["outputs/derived/upper_rr_no_CES/full*.tsv", "outputs/derived/upper_rr_no_CES/pooled_estimates.tsv", "outputs/derived/upper_rr_no_CES/UpperRR*"],
        "public_example": ["qa/public_example/*"], "korean": ["outputs/derived/korean_wgs/*", "outputs/figures/supplementary/SupplementaryFigure7_*"],
        "main": ["outputs/derived/main_autism/" + x for x in ("figure2_summary.rds", "mutational_variance.tsv", "figure2A_cohort_table.tsv", "fraction_cases_by_rr.tsv", "effective_penetrance.tsv", "gene_mutational_variance.tsv")],
        "cohort": ["outputs/derived/cohort_ddid/cohort_sex*"],
        "ddid": ["outputs/derived/cohort_ddid/ddid*"],
        "secondary": ["outputs/derived/main_autism/" + x for x in ("secondary_autism_summary.rds", "sibling_mutational_variance.tsv", "prevalence_mutational_variance.tsv")],
        "posterior": ["outputs/derived/main_autism/" + x for x in ("tada_comparison.tsv", "tada_comparison_statistics.tsv", "cds_by_rr.tsv", "posterior_diagnostics_summary.rds")],
        "prevalence": ["outputs/derived/autism_dd/prevalence_sweep.tsv"],
        "forecast": ["outputs/derived/forecasting/*"],
        "dd": ["outputs/derived/autism_dd/" + x for x in ("figure4_summary.rds", "mutational_variance.tsv", "fraction_cases.tsv", "effective_rr.tsv", "number_genes.tsv", "ddd_fraction_by_class.tsv")],
        "no_overlap": ["outputs/derived/sensitivities/no_overlap*"], "paternal": ["outputs/tables/supplementary/paternal_age/*", "outputs/figures/diagnostics/paternal_age/*"],
        "figure2": ["outputs/figures/main/Figure2.pdf"], "figure3": ["outputs/figures/main/Figure3.pdf", "outputs/figures/supplementary/SupplementaryFigure13_*"],
        "figure4": ["outputs/figures/main/Figure4.pdf"],
        "supplementary": [f"outputs/figures/supplementary/SupplementaryFigure{i}_*" for i in (2,3,4,5,6,8,9,10,11,12,14,15)],
        "tables": [f"outputs/tables/supplementary/SupplementaryTable{i}_*" for i in range(2,10)],
        "checklist": ["outputs/tables/manuscript_values/*"], "finalize": ["submission/**/*", "provenance/output_manifest.tsv"],
    }[stage]
    found = []
    for pattern in patterns:
        matching = [p for p in root.glob(pattern) if p.is_file()]
        if not matching:
            raise RuntimeError(f"Missing output for {stage}: {pattern}")
        found.extend(matching)
    found = sorted(set(found))
    if not found:
        raise RuntimeError(f"No outputs found for {stage}; check the output contract")
    return found


def record_path(c, stage, task=None):
    return Path(c["root"]) / "workflow/completed" / (stage + (f"_{task}" if task else "") + ".json")


def verify(c, stage, task=None):
    if stage in ARRAYS and task is None:
        for i in range(1, ARRAYS[stage] + 1):
            verify(c, stage, i)
        return
    record = json.loads(record_path(c, stage, task).read_text())
    if record["identity"] != sha(Path(c["root"]) / "workflow/identity.json"):
        raise RuntimeError(f"Stale stage: {stage}")
    for name, digest in record["outputs"].items():
        if sha(Path(c["root"]) / name) != digest:
            raise RuntimeError(f"Changed output: {name}")


def initialize(c):
    root = Path(c["root"])
    root.mkdir(parents=True, exist_ok=True)
    (root / "R-library").mkdir(exist_ok=True)
    identity = {"config": {k: v for k, v in c.items() if k not in ("slurm", "forecast_cores")}, "code": code_hashes()}
    path = root / "workflow/identity.json"
    if path.exists():
        if json.loads(path.read_text()) != identity:
            raise RuntimeError("Run configuration/code changed. Use a new root, not stale checkpoints.")
    else:
        if (root / "outputs").exists():
            raise RuntimeError("Unmanaged outputs already exist; use a new root")
        atomic_json(path, identity)
    atomic_json(root / "workflow/config.json", c)


def verify_inputs(c):
    """Before an explicit restart, reject changed raw inputs or saved draws."""
    root = Path(c["root"])
    for filename in ("input_manifest.tsv", "bootstrap_inputs.tsv"):
        path = root / "provenance" / filename
        if not path.exists():
            continue  # Imported reporting roots record imports instead of raw preparation.
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                source = Path(row["path"])
                if not source.is_absolute(): source = root / source
                if sha(source) != row["sha256"]:
                    raise RuntimeError("Changed input: " + str(source))


def execute(c, stage, task=None):
    root = Path(c["root"])
    initialize(c)
    if task is not None and (stage not in ARRAYS or not 1 <= task <= ARRAYS[stage]):
        raise ValueError("Invalid array task")
    for parent in PARENTS[stage]:
        verify(c, parent)
    destination = record_path(c, stage, task)
    lock = destination.with_suffix(".lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    lock.mkdir()  # Atomic, exclusive ownership; stale locks require inspection.
    try:
        if destination.exists():
            raise RuntimeError("Stage already completed; explicitly use --reuse")
        env = environment(c)
        if stage in ("prepare", "models"):
            env["BURDENMLEDN_JOB_STAGE"] = "prepare" if stage == "prepare" else "collect"
        cwd = root / "qa/public_example" if stage == "public_example" else REPO
        cwd.mkdir(parents=True, exist_ok=True)
        log = root / "logs" / (destination.stem + ".log")
        log.parent.mkdir(exist_ok=True)
        with log.open("w") as stream:
            for command in commands(stage, c, task):
                stream.write(shlex.join(command) + "\n"); stream.flush()
                subprocess.run(command, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
        outputs = {str(p.relative_to(root)): sha(p) for p in products(stage, c, task)}
        atomic_json(destination, {"identity": sha(root / "workflow/identity.json"), "outputs": outputs,
                                  "finished_unix": time.time()})
    finally:
        lock.rmdir()


def selection(names, stop_at=()):
    wanted = set()
    def add(name):
        if name not in PARENTS:
            raise ValueError(f"Unknown stage: {name}")
        if name not in stop_at:
            for parent in PARENTS[name]: add(parent)
        wanted.add(name)
    for name in names: add(name)
    return [name for name in PARENTS if name in wanted]


def submit(c, stages, reuse):
    root = Path(c["root"])
    ledger_path = root / "workflow/submitted.json"
    if ledger_path.exists():
        raise RuntimeError("Submission ledger exists; inspect previous jobs before another submission")
    jobs = {}
    for stage in stages:
        if stage in reuse:
            verify(c, stage); continue
        cpus = c["forecast_cores"] if stage == "forecast" else 1
        slurm = c["slurm"]
        args = ["sbatch", "--parsable", "--job-name=burden_" + stage,
                f"--cpus-per-task={cpus}", f"--mem={slurm['memory_gb']}G",
                "--partition=" + slurm["partition"], "--time=" + slurm["time"],
                "--output=" + str(root / f"logs/{stage}_%A_%a.out"),
                "--error=" + str(root / f"logs/{stage}_%A_%a.err")]
        parents = [jobs[p] for p in PARENTS[stage] if p in jobs]
        if parents: args += ["--dependency=afterok:" + ":".join(parents)]
        if stage in ARRAYS:
            limit = slurm["array_limit"]
            args += [f"--array=1-{ARRAYS[stage]}" + (f"%{int(limit)}" if limit else "")]
        command = [sys.executable, str(HERE / "run.py"), "--config", str(root / "workflow/config.json"), "--execute-stage", stage]
        if slurm["image"]:
            command = ["singularity", "exec", "--bind", f"{REPO}:{REPO}:ro,{root}:{root}:rw", *sum((["--bind", b] for b in slurm["binds"]), []), slurm["image"], "python3", *command[1:]]
        script = "#!/usr/bin/env bash\nset -euo pipefail\n" + shlex.join(command)
        if stage in ARRAYS: script += ' --task "$SLURM_ARRAY_TASK_ID"'
        result = subprocess.run(args, input=script + "\n", text=True, capture_output=True, check=True)
        job = result.stdout.strip().split(";")[0]
        if not job.isdigit(): raise RuntimeError("Unexpected sbatch response: " + result.stdout)
        jobs[stage] = job
        atomic_json(ledger_path, jobs)
        print(stage, job, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(HERE / "manuscript.json"))
    parser.add_argument("--backend", choices=("local", "slurm"), default="local")
    parser.add_argument("--stages", default="finalize", help="Comma-separated targets; includes dependencies")
    parser.add_argument("--reuse", default="", help="Comma-separated completed stages; checks hashes")
    parser.add_argument("--run", action="store_true", help="Execute/submit; otherwise print plan")
    parser.add_argument("--execute-stage", choices=PARENTS, help=argparse.SUPPRESS)
    parser.add_argument("--task", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    legacy = sorted(k for k in os.environ if k.startswith("RUN_") and k != "RUN_DATE")
    if legacy:
        parser.error("Legacy selectors are no longer supported: " + ", ".join(legacy) +
                     ". Unset them; use --stages and --reuse, inspect the plan, then add --run.")
    if len(sys.argv) == 1:
        parser.print_help(); return
    c = load_config(args.config)
    if args.execute_stage:
        execute(c, args.execute_stage, args.task); return
    reuse = set(filter(None, args.reuse.split(",")))
    stages = selection(args.stages.split(","), reuse)
    if not reuse <= set(stages): raise ValueError("Reuse stages must belong to selected graph")
    if not args.run:
        print("Stage\tDependencies\tTasks\tAction")
        for stage in stages:
            print(stage, ",".join(PARENTS[stage]) or "-", ARRAYS.get(stage, 1), "reuse" if stage in reuse else "run", sep="\t")
        return
    initialize(c)
    if reuse: verify_inputs(c)
    (Path(c["root"]) / "logs").mkdir(exist_ok=True)
    if args.backend == "slurm": submit(c, stages, reuse); return
    for stage in stages:
        if stage in reuse:
            verify(c, stage); continue
        for task in range(1, ARRAYS[stage] + 1) if stage in ARRAYS else [None]:
            execute(c, stage, task)


if __name__ == "__main__":
    main()
