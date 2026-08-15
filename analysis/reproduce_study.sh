#!/usr/bin/env bash
set -euo pipefail

# One-command entry point for reproducing the BurdenMLE-DN study.
# Run from any working directory with:
#   bash analysis/reproduce_study.sh
#
# Override Rscript if it is not on PATH:
#   R_BIN=/path/to/Rscript bash Final_Runs_July2026/repo/BurdenMLE_DN/reproduce_study.sh
#
# Stages are opt-in so an accidental invocation never starts a long model run.
# Example full figure rebuild from existing model outputs:
#   RUN_DATE=Aug15_26 RUN_SIMULATION_ANALYSIS=true \
#   RUN_MAIN_FIGURES=true RUN_FORECASTING=true \
#   RUN_AUTISM_DD_FIGURES=true RUN_SUPPLEMENTARY_FIGURES=true \
#   RUN_SUPPLEMENTARY_TABLES=true bash .../reproduce_study.sh

ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${ANALYSIS_DIR}/.." && pwd)"
SCRIPT_DIR="${ANALYSIS_DIR}/scripts"
FINAL_RUNS_DIR="${BURDENMLEDN_ANALYSIS_ROOT:-${ANALYSIS_DIR}}"
LOG_DIR="${FINAL_RUNS_DIR}/outputs/logs"
R_BIN="${R_BIN:-Rscript}"
export BURDENMLEDN_ANALYSIS_ROOT="${FINAL_RUNS_DIR}"

if ! command -v "${R_BIN}" >/dev/null 2>&1; then
  MACOS_RSCRIPT="/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/Rscript"
  if [[ -x "${MACOS_RSCRIPT}" ]]; then
    R_BIN="${MACOS_RSCRIPT}"
  else
    echo "Could not find Rscript. Set R_BIN to its absolute path." >&2
    exit 1
  fi
fi

RUN_DATE="${RUN_DATE:-$(date '+%b%d_%y')}"
MAIN_MODEL_MANIFEST="${MAIN_MODEL_MANIFEST:-${FINAL_RUNS_DIR}/outputs/data/model_manifest_main_${RUN_DATE}.rds}"
NO_CES_MODEL_MANIFEST="${NO_CES_MODEL_MANIFEST:-${FINAL_RUNS_DIR}/outputs/data/model_manifest_no_ces_${RUN_DATE}.rds}"
NO_OVERLAP_MODEL_MANIFEST="${NO_OVERLAP_MODEL_MANIFEST:-${FINAL_RUNS_DIR}/outputs/data/model_manifest_no_overlap_${RUN_DATE}.rds}"
SIMULATION_SEED="${SIMULATION_SEED:-20260715}"
BOOTSTRAPS="${BOOTSTRAPS:-100}"
MODEL_SEED="${MODEL_SEED:-24312342}"
KOREAN_WGS_SEED="${KOREAN_WGS_SEED:-20260721}"
KOREAN_WGS_BOOTSTRAPS="${KOREAN_WGS_BOOTSTRAPS:-${BOOTSTRAPS}}"
if [[ ! "${RUN_DATE}" =~ ^[A-Za-z0-9_-]+$ ]]; then
  echo "RUN_DATE may contain only letters, numbers, underscores, and hyphens." >&2
  exit 1
fi
# Every stage defaults to false; select only the layers that need rebuilding.
RUN_ALL="${RUN_ALL:-false}"
RUN_SIMULATIONS="${RUN_SIMULATIONS:-${RUN_ALL}}"
RUN_EM_SIMULATIONS="${RUN_EM_SIMULATIONS:-${RUN_SIMULATIONS}}"
RUN_MIXSQP_SIMULATIONS="${RUN_MIXSQP_SIMULATIONS:-${RUN_SIMULATIONS}}"
RUN_SIMULATION_ANALYSIS="${RUN_SIMULATION_ANALYSIS:-${RUN_ALL}}"
RUN_KOREAN_WGS="${RUN_KOREAN_WGS:-${RUN_ALL}}"
RUN_MAIN_MODELS="${RUN_MAIN_MODELS:-${RUN_ALL}}"
RUN_NO_CES_MODELS="${RUN_NO_CES_MODELS:-${RUN_ALL}}"
RUN_NO_OVERLAP_MODELS="${RUN_NO_OVERLAP_MODELS:-${RUN_ALL}}"
RUN_MAIN_FIGURES="${RUN_MAIN_FIGURES:-${RUN_ALL}}"
RUN_FORECASTING="${RUN_FORECASTING:-${RUN_ALL}}"
RUN_AUTISM_DD_FIGURES="${RUN_AUTISM_DD_FIGURES:-${RUN_ALL}}"
RUN_SUPPLEMENTARY_FIGURES="${RUN_SUPPLEMENTARY_FIGURES:-${RUN_ALL}}"
RUN_SUPPLEMENTARY_TABLES="${RUN_SUPPLEMENTARY_TABLES:-${RUN_ALL}}"
RUN_MANUSCRIPT_ESTIMATES="${RUN_MANUSCRIPT_ESTIMATES:-${RUN_ALL}}"

for stage_flag in RUN_ALL RUN_EM_SIMULATIONS RUN_MIXSQP_SIMULATIONS RUN_SIMULATION_ANALYSIS RUN_KOREAN_WGS RUN_MAIN_MODELS RUN_NO_CES_MODELS RUN_NO_OVERLAP_MODELS RUN_MAIN_FIGURES RUN_FORECASTING RUN_AUTISM_DD_FIGURES RUN_SUPPLEMENTARY_FIGURES RUN_SUPPLEMENTARY_TABLES RUN_MANUSCRIPT_ESTIMATES; do
  if [[ "${!stage_flag}" != "true" && "${!stage_flag}" != "false" ]]; then
    echo "${stage_flag} must be true or false." >&2
    exit 1
  fi
done

mkdir -p "${LOG_DIR}"

run_step() {
  local step_name="$1"
  shift
  local log_file="${LOG_DIR}/${RUN_DATE}_${step_name}.log"

  echo
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting ${step_name}"
  "$@" 2>&1 | tee "${log_file}"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished ${step_name}"
}

if [[ "${RUN_EM_SIMULATIONS}" == "true" ]]; then
  run_step "simulation_em_seed${SIMULATION_SEED}" \
    "${R_BIN}" "${SCRIPT_DIR}/simulation_study.R" \
    --optimizer EM \
    --run-date "${RUN_DATE}" \
    --seed "${SIMULATION_SEED}"

fi

if [[ "${RUN_MIXSQP_SIMULATIONS}" == "true" ]]; then
  run_step "simulation_mixsqp_seed${SIMULATION_SEED}" \
    "${R_BIN}" "${SCRIPT_DIR}/simulation_study.R" \
    --optimizer mixsqp \
    --run-date "${RUN_DATE}" \
    --seed "${SIMULATION_SEED}"
fi

if [[ "${RUN_SIMULATION_ANALYSIS}" == "true" ]]; then
  run_step "analyze_simulations_mixsqp_vs_em_seed${SIMULATION_SEED}" \
    "${R_BIN}" "${SCRIPT_DIR}/analyze_simulations.R" \
    --optimizer mixsqp \
    --compare-optimizer EM \
    --run-date "${RUN_DATE}" \
    --seed "${SIMULATION_SEED}"
fi

if [[ "${RUN_KOREAN_WGS}" == "true" ]]; then
  run_step "run_korean_wgs_mixsqp" \
    "${R_BIN}" "${SCRIPT_DIR}/run_korean_wgs.R" \
    --run-date "${RUN_DATE}" \
    --seed "${KOREAN_WGS_SEED}" \
    --n-boot "${KOREAN_WGS_BOOTSTRAPS}"
fi

if [[ "${RUN_MAIN_MODELS}" == "true" ]]; then
  run_step "run_models_main_mixsqp" \
    "${R_BIN}" "${SCRIPT_DIR}/run_models.R" \
    --mode main \
    --optimizer mixsqp \
    --run-autism true \
    --run-ddd true \
    --run-date "${RUN_DATE}" \
    --seed "${MODEL_SEED}" \
    --bootstraps "${BOOTSTRAPS}"
fi

if [[ "${RUN_NO_CES_MODELS}" == "true" ]]; then
  run_step "run_models_no_ces_mixsqp" \
    "${R_BIN}" "${SCRIPT_DIR}/run_models.R" \
    --mode no_ces \
    --optimizer mixsqp \
    --run-autism true \
    --run-ddd false \
    --run-date "${RUN_DATE}" \
    --seed "${MODEL_SEED}" \
    --bootstraps "${BOOTSTRAPS}"
fi

if [[ "${RUN_NO_OVERLAP_MODELS}" == "true" ]]; then
  run_step "run_models_no_overlap_mixsqp" \
    "${R_BIN}" "${SCRIPT_DIR}/run_models.R" \
    --mode no_overlap \
    --optimizer mixsqp \
    --run-autism true \
    --run-ddd false \
    --run-date "${RUN_DATE}" \
    --seed "${MODEL_SEED}" \
    --bootstraps "${BOOTSTRAPS}"
fi

if [[ "${RUN_MAIN_FIGURES}" == "true" ]]; then
  run_step "summarize_main_autism" \
    "${R_BIN}" "${SCRIPT_DIR}/summarize_main_autism.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
  run_step "summarize_cohort_sex" \
    "${R_BIN}" "${SCRIPT_DIR}/summarize_cohort_sex.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
  run_step "make_figure2" \
    "${R_BIN}" "${SCRIPT_DIR}/make_figure2_cohort_sex.R"
fi

if [[ "${RUN_FORECASTING}" == "true" ]]; then
  run_step "forecast_gene_discovery" \
    "${R_BIN}" "${SCRIPT_DIR}/forecasting_script_revision.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}" \
    --seed "${MODEL_SEED}" \
    --cores "${FORECAST_CORES:-1}"
  run_step "make_forecasting_figures" \
    "${R_BIN}" "${SCRIPT_DIR}/make_forecasting_figures.R"
fi

if [[ "${RUN_AUTISM_DD_FIGURES}" == "true" ]]; then
  run_step "summarize_autism_dd" \
    "${R_BIN}" "${SCRIPT_DIR}/summarize_autism_dd.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
  run_step "make_figure4_autism_dd" \
    "${R_BIN}" "${SCRIPT_DIR}/make_figure4_autism_dd.R"
  run_step "summarize_autism_dd_prevalence" \
    "${R_BIN}" "${SCRIPT_DIR}/summarize_autism_dd_prevalence.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
fi

if [[ "${RUN_SUPPLEMENTARY_FIGURES}" == "true" ]]; then
  run_step "summarize_secondary_autism" \
    "${R_BIN}" "${SCRIPT_DIR}/summarize_secondary_autism.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
  run_step "summarize_posterior_diagnostics" \
    "${R_BIN}" "${SCRIPT_DIR}/summarize_posterior_diagnostics.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
  run_step "summarize_ddid" \
    "${R_BIN}" "${SCRIPT_DIR}/summarize_ddid.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}" \
    --cores "${SECONDARY_CORES:-1}"
  run_step "summarize_no_overlap" \
    "${R_BIN}" "${SCRIPT_DIR}/summarize_no_overlap.R" \
    --model-manifest "${NO_OVERLAP_MODEL_MANIFEST}"
  run_step "describe_us_father_age" \
    "${R_BIN}" "${SCRIPT_DIR}/describe_us_father_age.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
  run_step "compare_spark_us_father_age" \
    "${R_BIN}" "${SCRIPT_DIR}/compare_spark_us_father_age.R"
  run_step "make_supplementary_figures" \
    "${R_BIN}" "${SCRIPT_DIR}/make_supplementary_figures.R" \
    --no-ces-model-manifest "${NO_CES_MODEL_MANIFEST}"
fi

if [[ "${RUN_SUPPLEMENTARY_TABLES}" == "true" ]]; then
  run_step "make_supplementary_tables" \
    "${R_BIN}" "${SCRIPT_DIR}/make_supplementary_tables.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
fi

if [[ "${RUN_MANUSCRIPT_ESTIMATES}" == "true" ]]; then
  run_step "make_manuscript_estimate_table" \
    "${R_BIN}" "${SCRIPT_DIR}/make_manuscript_estimate_table.R" \
    --model-manifest "${MAIN_MODEL_MANIFEST}"
fi

echo
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Requested study stages completed for ${RUN_DATE}."
