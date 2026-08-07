#!/usr/bin/env bash
# =============================================================================
# CrI3-pbc-ucf -> mumaxplus Moire workflow
#
# mumaxross Python env:  conda activate mumaxross
# GPU memcheck (optional): ./run_mumax_memcheck.sh …
#   Default CUDA tool: COMPUTE_SANITIZER=/usr/local/cuda-12.8/compute-sanitizer/compute-sanitizer
# =============================================================================
# C++ ./main writes twisted-double-bilayer coarse output by default
# (moire_twisted_double_bilayer_coarse.bin). Legacy twisted-bilayer
# (moire_coarse_v2.bin) is omitted unless you pass ``--moire-export both`` to ./main.
#
# This script runs the twisted-double-bilayer pipeline by default; set
# RUN_TWISTED_BILAYER=1 if you also need the legacy NPZ/plots from moire_coarse_v2.bin.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# --- configurable paths (override on command line: VAR=value ./run_moire_to_mumaxplus.sh)
: "${RUN_MAIN:=0}"
: "${DO_PLOT:=1}"
: "${DO_MUMAX:=1}"
: "${MUMAX_VERIFY_PLOT:=1}"
: "${PERIODS_X:=2}"
: "${PERIODS_Y:=4}"
# Twisted-double-bilayer: MUMAX_EXCHANGE_ONLY=1 → crI3_twisted_double_bilayer_moire_mumaxplus_setup.py --exchange-only
: "${MUMAX_EXCHANGE_ONLY:=1}"
# Twisted-double-bilayer: DMI_WORKFLOW_PLOT=1 → plot_twisted_double_bilayer_moire_exchange_npz.py --dmi-workflow
: "${DMI_WORKFLOW_PLOT:=0}"

: "${RUN_TWISTED_BILAYER:=0}"
: "${RUN_TWISTED_DOUBLE_BILAYER:=1}"

# Preserve the original variable names for the active twisted-bilayer path.
: "${MOIRE_BIN:=moire_coarse_v2.bin}"
: "${MOIRE_NPZ:=twisted_bilayer_moire_exchange.npz}"
: "${FIG_PREFIX:=twisted_bilayer_moire_workflow}"

: "${TWISTED_DOUBLE_BILAYER_BIN:=moire_twisted_double_bilayer_coarse.bin}"
: "${TWISTED_DOUBLE_BILAYER_NPZ:=twisted_double_bilayer_moire_exchange.npz}"
: "${TWISTED_DOUBLE_BILAYER_FIG_PREFIX:=twisted_double_bilayer_moire_workflow}"

# --- Example ./main invocation (twist angle [deg], max inter range [A], ...)
MAIN_EXAMPLE=(1.1 9.99)
# Optional extra args: "" inter_scale dmi_inter_scale j_twist j_prist j_intra dmi_sub dmi_basis
# MAIN_EXAMPLE=(1.1 9.99 "" 1.0 1.0 1.0 1.0 1.0 1.0 3)

run_pipeline() {
  local label="$1"
  local bin_path="$2"
  local npz_path="$3"
  local fig_prefix="$4"
  local bin_to_npz_script="$5"
  local validate_script="$6"
  local plot_script="$7"
  local mumax_script="$8"

  if [[ ! -f "${bin_path}" ]]; then
    echo "Missing ${bin_path} for ${label} (run ./main first, or set ${label} path vars)." >&2
    return 1
  fi

  echo "==> [${label}] Converting ${bin_path} -> ${npz_path}"
  python3 "${bin_to_npz_script}" "${bin_path}" -o "${npz_path}"

  echo "==> [${label}] Validating ${npz_path}"
  python3 "${validate_script}" --validate "${npz_path}"

  if [[ "${DO_PLOT}" == "1" ]]; then
    echo "==> [${label}] Plotting maps -> ${fig_prefix}_*.png"
    export MPLBACKEND="${MPLBACKEND:-Agg}"
    if [[ "${label}" == "twisted-bilayer" ]]; then
      python3 "${plot_script}" \
        --verify-npz "${npz_path}" \
        --show-converted \
        --save-prefix "${fig_prefix}"
    elif [[ "${label}" == "twisted-double-bilayer" ]]; then
      local -a plot_args=(
        python3 "${plot_script}"
        --verify-npz "${npz_path}"
        --save-prefix "${fig_prefix}"
        --periods-x "${PERIODS_X}"
        --periods-y "${PERIODS_Y}"
      )
      if [[ "${MUMAX_EXCHANGE_ONLY}" == "1" ]]; then
        plot_args+=(--exchange-only)
      fi
      if [[ "${DMI_WORKFLOW_PLOT}" == "1" ]]; then
        plot_args+=(--dmi-workflow)
      fi
      "${plot_args[@]}"
    else
      python3 "${plot_script}" \
        --verify-npz "${npz_path}" \
        --save-prefix "${fig_prefix}"
    fi
  fi

  if [[ "${DO_MUMAX}" == "1" ]]; then
    local -a mumax_args=(
      python3 "${mumax_script}" "${npz_path}"
      --save-prefix "${fig_prefix}_mumax"
      --periods-x "${PERIODS_X}"
      --periods-y "${PERIODS_Y}"
    )
    if [[ "${MUMAX_VERIFY_PLOT}" == "1" ]]; then
      mumax_args+=(--verify-plot)
    fi
    if [[ "${label}" == "twisted-double-bilayer" && "${MUMAX_EXCHANGE_ONLY}" == "1" ]]; then
      mumax_args+=(--exchange-only)
    fi
    echo "==> [${label}] mumaxplus setup: ${mumax_args[*]}"
    "${mumax_args[@]}"
  fi
}

# -----------------------------------------------------------------------------
# Step 0 (optional): run atomistic binary
# -----------------------------------------------------------------------------
if [[ "${RUN_MAIN}" == "1" ]]; then
  if [[ ! -x ./main ]]; then
    echo "RUN_MAIN=1 but ./main missing; run: make -j4" >&2
    exit 1
  fi
  echo "Running ./main ${MAIN_EXAMPLE[*]} ..."
  ./main "${MAIN_EXAMPLE[@]}"
fi

if [[ "${RUN_TWISTED_BILAYER}" == "1" ]]; then
  run_pipeline \
    "twisted-bilayer" \
    "${MOIRE_BIN}" \
    "${MOIRE_NPZ}" \
    "${FIG_PREFIX}" \
    "twisted_bilayer_moire_coarse_bin_to_npz.py" \
    "twisted_bilayer_moire_exchange_convert.py" \
    "plot_twisted_bilayer_moire_exchange_npz.py" \
    "crI3_twisted_bilayer_moire_mumaxplus_setup.py"
fi

if [[ "${RUN_TWISTED_DOUBLE_BILAYER}" == "1" ]]; then
  run_pipeline \
    "twisted-double-bilayer" \
    "${TWISTED_DOUBLE_BILAYER_BIN}" \
    "${TWISTED_DOUBLE_BILAYER_NPZ}" \
    "${TWISTED_DOUBLE_BILAYER_FIG_PREFIX}" \
    "twisted_double_bilayer_moire_coarse_bin_to_npz.py" \
    "twisted_double_bilayer_moire_exchange_convert.py" \
    "plot_twisted_double_bilayer_moire_exchange_npz.py" \
    "crI3_twisted_double_bilayer_moire_mumaxplus_setup.py"
fi

echo "Done. Key outputs:"
if [[ "${RUN_TWISTED_BILAYER}" == "1" ]]; then
  echo "  ${MOIRE_NPZ} - twisted-bilayer micromagnetic input"
  echo "  ${FIG_PREFIX}_*.png - twisted-bilayer plots"
fi
if [[ "${RUN_TWISTED_DOUBLE_BILAYER}" == "1" ]]; then
  echo "  ${TWISTED_DOUBLE_BILAYER_NPZ} - twisted-double-bilayer micromagnetic input"
  echo "  ${TWISTED_DOUBLE_BILAYER_FIG_PREFIX}_*.png - twisted-double-bilayer plots"
fi
