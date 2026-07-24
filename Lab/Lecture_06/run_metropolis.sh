#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# NSL Exercise 06.1 - Automatic temperature scan (Metropolis)
# T: 2.0 -> 0.5 step -0.05
# ============================================================

# Directory in cui si trova questo script 
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SIM_DIR="$BASE_DIR/ex06.1"
SRC_DIR="$SIM_DIR/SOURCE"
IN_DIR="$SIM_DIR/INPUT"
OUT_DIR="$SIM_DIR/OUTPUT"

INPUT_DAT="$IN_DIR/input.dat"
PROPS_DAT="$IN_DIR/properties.dat"
ARCHIVE_BASE="$SIM_DIR/OUTPUT_metropolis"

# Parametri simulazione
SIM_TYPE=2          # 2 = Metropolis, 3 = Gibbs
JCOUP=1.0
NPART=50
RHO=1.0
RCUT=0
DELTA=0
NBLOCKS=20
NSTEPS=20000

FIRST_RUN=1

echo ">>> Compiling..."
make -C "$SRC_DIR" clean
make -C "$SRC_DIR"

mkdir -p "$ARCHIVE_BASE"

write_input_dat () {
  local temp="$1"
  local restart="$2"
  local hfield="$3"

  cat > "$INPUT_DAT" <<EOF
SIMULATION_TYPE        ${SIM_TYPE}    ${JCOUP}    ${hfield}
RESTART                ${restart}
TEMP                   ${temp}
NPART                  ${NPART}
RHO                    ${RHO}
R_CUT                  ${RCUT}
DELTA                  ${DELTA}
NBLOCKS                ${NBLOCKS}
NSTEPS                 ${NSTEPS}

ENDINPUT
EOF
}

write_properties_h0 () {
  cat > "$PROPS_DAT" <<EOF
TOTAL_ENERGY
SPECIFIC_HEAT
SUSCEPTIBILITY

ENDPROPERTIES
EOF
}

write_properties_h002 () {
  cat > "$PROPS_DAT" <<EOF
MAGNETIZATION

ENDPROPERTIES
EOF
}

run_and_archive () {
  local temp="$1"
  local hlabel="$2"
  local obslabel="$3"

  local tdir="$ARCHIVE_BASE/T_${temp}"
  local rdir="$tdir/${hlabel}_${obslabel}"
  mkdir -p "$rdir"

  echo ">>> Run T=${temp}, ${hlabel}, ${obslabel}"
  ( cd "$SRC_DIR" && ./simulator.exe )

  cp -r "$OUT_DIR"/* "$rdir"/
  cp "$INPUT_DAT" "$rdir/input_used.dat"
  cp "$PROPS_DAT" "$rdir/properties_used.dat"
}

for T in $(awk 'BEGIN{for(t=2.00;t>=0.50-1e-9;t-=0.05) printf "%.2f\n", t}'); do

  # Run 1: h=0.00 per U, C, chi
  if [[ "$FIRST_RUN" -eq 1 ]]; then
    RESTART=0
  else
    RESTART=1
  fi

  write_input_dat "$T" "$RESTART" "0.00"
  write_properties_h0
  run_and_archive "$T" "h0.00" "thermodynamics"

  FIRST_RUN=0

  # Run 2: h=0.02 per M
  write_input_dat "$T" "1" "0.02"
  write_properties_h002
  run_and_archive "$T" "h0.02" "magnetization"

done

echo ">>> Scan completed."
echo "Risults in: $ARCHIVE_BASE"

# per renderlo eseguibile: $ chmod +x run_metropolis.sh
# per avviarlo: ./run_metropolis.sh