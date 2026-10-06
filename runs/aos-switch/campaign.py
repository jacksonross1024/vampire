#!/usr/bin/env python3
"""Run aos-switch cases, at most two MPI jobs at a time.

Each case is a copy of input-aos-sv and CoPt-aos-sv.mat. The normal edits
are the Co-1 and Co-2 initial spin, sim:laser-pulse-time,
sim:laser-pulse-power, and sim:integrator-random-seed. Interface
conductances are rewritten to the finite values below. 1e12 m/s is the
transparent comparison limit and is not the default. Setting
AOS_DISABLE_XY_PBC=1 also strips the x/y periodic boundary lines from the
per-run input copy for finite-lateral-size controls.
Co-1 and Co-2 are the thin cobalt above the Cu. Their initial spin is
0,0,+1 or 0,0,-1. Co-3..Co-6 are the thick reference under the Cu and
stay at 0,0,1. Fluence and pulse duration are the other scan axes.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = Path(__file__).resolve().parent
RUN_ROOT = Path(os.environ.get("AOS_RUN_ROOT", CAMPAIGN / "newstack"))
BINARY = ROOT / "vampire-parallel"
INPUT = ROOT / "input-aos-sv"
MAT = ROOT / "CoPt-aos-sv.mat"
ABS = ROOT / "abs.txt"
DISABLE_XY_PBC = os.environ.get("AOS_DISABLE_XY_PBC", "").lower() in {"1", "true", "yes"}
CPU_SETS = ["0-7", "8-15"]
MAX_JOBS = 2

# Co-1, Co-2: materials 3 and 5. Co-3..Co-6: materials 7, 9, 11, 13.
TOP_CO = (3, 5)
REF_CO = (7, 9, 11, 13)

# Finite interface conductance (m/s): material, partner above it, G.
# Co/Pt and Ta/Pt 2e4, Co/Cu 3e4, Ta/SiO2 0.1. Same block as CoPt-aos-sv.mat.
INTERFACE_CONDUCTANCE = (
    (16, 15, "0.1", "SiO2/Ta"),
    (15, 14, "2.0e4", "Ta/Pt"),
    (14, 13, "2.0e4", "Pt/Co"),
    (13, 12, "2.0e4", "Co/Pt"),
    (12, 11, "2.0e4", "Pt/Co"),
    (11, 10, "2.0e4", "Co/Pt"),
    (10, 9, "2.0e4", "Pt/Co"),
    (9, 8, "2.0e4", "Co/Pt"),
    (8, 7, "2.0e4", "Pt/Co"),
    (7, 6, "3.0e4", "Co/Cu"),
    (6, 5, "3.0e4", "Cu/Co"),
    (5, 4, "2.0e4", "Co/Pt"),
    (4, 3, "2.0e4", "Pt/Co"),
    (3, 2, "2.0e4", "Co/Pt"),
    (2, 1, "2.0e4", "Pt/Ta"),
)


def tag_pw(picoseconds):
    return f"pw{picoseconds:.4g}ps"


def tag_fluence(fluence):
    return f"F{fluence:.4g}"


def tag_mz(top_mz):
    return "mz-1" if top_mz < 0 else "mz+1"


def case_dir(pulse_ps, fluence, top_mz, seed):
    return RUN_ROOT / tag_pw(pulse_ps) / tag_fluence(fluence) / tag_mz(top_mz) / f"seed{seed:02d}"


def case_id(pulse_ps, fluence, top_mz, seed):
    return f"{tag_pw(pulse_ps)}/{tag_fluence(fluence)}/{tag_mz(top_mz)}/seed{seed:02d}"


def apply_interface_conductance(mat):
    """Write the finite interface conductances onto the copied mat."""
    for material, partner, conductance, note in INTERFACE_CONDUCTANCE:
        mat, nsub = re.subn(
            rf"^material\[{material}\]:spin-interface-conductance\s*=.*$",
            f"material[{material}]:spin-interface-conductance = {partner}, {conductance}       # {note}",
            mat,
            count=1,
            flags=re.M,
        )
        if nsub != 1:
            raise RuntimeError(f"failed to set interface conductance for material {material}")
    return mat


def write_case(pulse_ps, fluence, top_mz, seed):
    folder = case_dir(pulse_ps, fluence, top_mz, seed)
    folder.mkdir(parents=True, exist_ok=True)
    text = INPUT.read_text()
    if DISABLE_XY_PBC:
        text = re.sub(r"^create:periodic-boundaries-[xy]\s*\n?", "", text, flags=re.M)
    text = re.sub(
        r"^sim:laser-pulse-power\s*=.*$",
        f"sim:laser-pulse-power = {fluence:.8g}",
        text,
        count=1,
        flags=re.M,
    )
    text = re.sub(
        r"^sim:laser-pulse-time\s*=.*$",
        f"sim:laser-pulse-time = {pulse_ps:.8g} !ps",
        text,
        count=1,
        flags=re.M,
    )
    text = re.sub(r"^sim:integrator-random-seed\s*=.*\n?", "", text, flags=re.M)
    text = text.rstrip() + f"\nsim:integrator-random-seed = {seed}\n"
    topology = os.environ.get("AOS_MPI_TOPOLOGY", "").strip()
    if topology:
        text, nsub = re.subn(
            r"^sim:mpi-topology\s*=.*$",
            f"sim:mpi-topology = {topology}",
            text,
            count=1,
            flags=re.M,
        )
        if nsub != 1:
            raise RuntimeError("failed to set sim:mpi-topology")
    (folder / "input").write_text(text)

    mat = MAT.read_text()
    direction = f"0,0,{int(top_mz)}"
    for index in TOP_CO:
        mat, nsub = re.subn(
            rf"^(material\[{index}\]:initial-spin-direction\s*=).*$",
            rf"\1 {direction}",
            mat,
            count=1,
            flags=re.M,
        )
        if nsub != 1:
            raise RuntimeError(f"failed to set Co material {index} spin")
    for index in REF_CO:
        mat, nsub = re.subn(
            rf"^(material\[{index}\]:initial-spin-direction\s*=).*$",
            r"\1 0,0,1",
            mat,
            count=1,
            flags=re.M,
        )
        if nsub != 1:
            raise RuntimeError(f"failed to set Co material {index} spin")
    mat = apply_interface_conductance(mat)
    (folder / "CoPt-aos-sv.mat").write_text(mat)
    shutil.copy(ABS, folder / "abs.txt")
    # ofstream does not create this directory, so the spin profiles are lost without it.
    (folder / "spin-acc").mkdir(exist_ok=True)
    meta = {
        "id": case_id(pulse_ps, fluence, top_mz, seed),
        "path": str(folder),
        "pulse_duration_ps": pulse_ps,
        "fluence_J_m2": fluence,
        "co12_mz": int(top_mz),
        "seed": int(seed),
        "xy_periodic_boundaries": not DISABLE_XY_PBC,
    }
    (folder / "case.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def finished(meta):
    log = Path(meta["path"]) / "run.log"
    if not log.exists():
        return False
    text = log.read_text(errors="replace")
    return "Simulation ended gracefully" in text


def mpi_ranks():
    return os.environ.get("AOS_NP", "8")


def cpu_sets():
    raw = os.environ.get("AOS_CPU_SETS", "").strip()
    if not raw:
        return list(CPU_SETS)
    return [part.strip() for part in raw.split(";") if part.strip()]


def launch(meta, cpu_set):
    folder = Path(meta["path"])
    log = open(folder / "run.log", "w")
    command = ["mpirun", "-np", mpi_ranks()]
    bind = os.environ.get("AOS_MPI_BIND", "1").lower() not in {"0", "false", "no"}
    if bind:
        command += ["--bind-to", "core", "--cpu-set", cpu_set]
    command.append(str(BINARY))
    print(f"START {meta['id']} cpu={cpu_set}", flush=True)
    process = subprocess.Popen(
        command,
        cwd=folder,
        stdout=log,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    return process, log, time.time()


def prune_profiles(meta):
    """Keep magnetisation output. Spin-current and temperature files repeat across seeds."""
    if int(meta["seed"]) == 1:
        return
    folder = Path(meta["path"])
    spin_dir = folder / "spin-acc"
    if spin_dir.is_dir():
        shutil.rmtree(spin_dir)
    temperature = folder / "vertical_temperature_profile.dat"
    if temperature.exists():
        temperature.unlink()


def analyse(meta, wall_time_s):
    payload = dict(meta)
    payload["wall_time_s"] = wall_time_s
    subprocess.check_call(
        [sys.executable, str(CAMPAIGN / "analyse.py"), json.dumps(payload)],
        cwd=str(CAMPAIGN),
    )
    prune_profiles(meta)


def run_queue(jobs):
    pending = [job for job in jobs if not finished(job)]
    active = []
    used_cpus = set()
    while pending or active:
        while pending and len(active) < MAX_JOBS:
            cpu = next(cpu for cpu in cpu_sets() if cpu not in used_cpus)
            meta = pending.pop(0)
            process, log, started = launch(meta, cpu)
            used_cpus.add(cpu)
            active.append((process, log, started, meta, cpu))
        if not active:
            break
        time.sleep(10)
        still = []
        for process, log, started, meta, cpu in active:
            code = process.poll()
            if code is None:
                still.append((process, log, started, meta, cpu))
                continue
            log.close()
            used_cpus.discard(cpu)
            wall = time.time() - started
            if code != 0 or not finished(meta):
                print(f"FAILED {meta['id']} exit={code}", flush=True)
                continue
            print(f"FINISH {meta['id']} wall_s={wall:.0f}", flush=True)
            analyse(meta, wall)
        active = still
        free_gb = shutil.disk_usage(CAMPAIGN).free / 1e9
        if free_gb < 8.0:
            print(f"STOP low disk free_gb={free_gb:.1f}", flush=True)
            for process, log, _, meta, _ in active:
                process.terminate()
                log.close()
            raise SystemExit(2)
    print("QUEUE_EMPTY", flush=True)


def baseline_jobs():
    jobs = []
    for top_mz in (-1, 1):
        for seed in range(1, 11):
            jobs.append(write_case(0.1, 55.0, top_mz, seed))
    return jobs


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "baseline"
    if mode == "baseline":
        run_queue(baseline_jobs())
        print("CAMPAIGN_BASELINE_COMPLETE", flush=True)
    elif mode == "queue":
        specs = json.loads(Path(sys.argv[2]).read_text())
        jobs = [
            write_case(
                item["pulse_ps"],
                item["fluence"],
                item["co12_mz"] if "co12_mz" in item else item["co14_mz"],
                item["seed"],
            )
            for item in specs
        ]
        run_queue(jobs)
        print("CAMPAIGN_QUEUE_COMPLETE", flush=True)
    else:
        raise SystemExit(f"unknown mode {mode}")
