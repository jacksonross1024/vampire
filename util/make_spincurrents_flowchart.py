#!/usr/bin/env python3
"""
Generate a diagrammatic flowchart (PNG) for the algorithm in:
  src/spintorque/spincurrents.cpp

No external dependencies beyond matplotlib.

Output:
  spincurrents-flowchart.png (in repo root by default)
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import List, Tuple, Optional

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


@dataclass
class Box:
    id: str
    text: str
    xy: Tuple[float, float]  # center (x,y) in figure coords
    wh: Tuple[float, float]  # width, height
    fc: str = "#F8FAFC"
    ec: str = "#0F172A"
    lw: float = 1.4


@dataclass
class Edge:
    src: str
    dst: str
    label: Optional[str] = None
    color: str = "#0F172A"


def draw_box(ax, b: Box):
    x, y = b.xy
    w, h = b.wh
    rect = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=b.lw,
        edgecolor=b.ec,
        facecolor=b.fc,
    )
    ax.add_patch(rect)
    ax.text(
        x,
        y,
        b.text,
        ha="center",
        va="center",
        fontsize=9.5,
        family="DejaVu Sans",
        color="#0B1220",
        wrap=True,
    )


def _box_anchor(box: Box, which: str) -> Tuple[float, float]:
    x, y = box.xy
    w, h = box.wh
    if which == "top":
        return (x, y + h / 2)
    if which == "bottom":
        return (x, y - h / 2)
    if which == "left":
        return (x - w / 2, y)
    if which == "right":
        return (x + w / 2, y)
    raise ValueError(which)


def draw_edge(ax, boxes_by_id, e: Edge):
    s = boxes_by_id[e.src]
    d = boxes_by_id[e.dst]

    # Simple routing heuristic based on relative positions.
    if d.xy[1] < s.xy[1]:
        p1 = _box_anchor(s, "bottom")
        p2 = _box_anchor(d, "top")
    elif d.xy[1] > s.xy[1]:
        p1 = _box_anchor(s, "top")
        p2 = _box_anchor(d, "bottom")
    elif d.xy[0] > s.xy[0]:
        p1 = _box_anchor(s, "right")
        p2 = _box_anchor(d, "left")
    else:
        p1 = _box_anchor(s, "left")
        p2 = _box_anchor(d, "right")

    arr = FancyArrowPatch(
        p1,
        p2,
        arrowstyle="-|>",
        mutation_scale=13,
        linewidth=1.2,
        color=e.color,
        connectionstyle="arc3,rad=0.0",
    )
    ax.add_patch(arr)

    if e.label:
        lx = (p1[0] + p2[0]) / 2
        ly = (p1[1] + p2[1]) / 2
        ax.text(
            lx,
            ly,
            e.label,
            ha="center",
            va="center",
            fontsize=8.5,
            color="#334155",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85),
        )


def build_flow() -> Tuple[List[Box], List[Edge]]:
    """
    High-level flow extracted from spincurrents.cpp:
    - initialise_spincurrents_1d(): allocate/size arrays per local MPI stacks
    - calculate_spin_accumulation_lepadatu_1d(): per time step, per local stack:
        charge transient (ns/ne) implicit step + compute Jc
        spin accumulation S(z,t) via Strang splitting:
            diffusion(dt/2) implicit
            explicit Heun for source/relax/precession
            diffusion(dt/2) implicit
        coarse average to microcells; compute torque field; fill outputs
      then MPI Allreduce for global arrays
    """

    boxes = [
        Box(
            id="init",
            text="initialise_spincurrents_1d()\nAllocate fine-grid arrays\nper local (MPI) stack\n(set default m̂_ref, etc.)",
            xy=(0.20, 0.88),
            wh=(0.34, 0.14),
            fc="#EEF2FF",
        ),
        Box(
            id="entry",
            text="calculate_spin_accumulation_lepadatu_1d()\n(per LLG step)",
            xy=(0.50, 0.88),
            wh=(0.40, 0.11),
            fc="#ECFEFF",
        ),
        Box(
            id="clear",
            text="Clear global outputs\nspin_torque, sa_final\nns/ne/jc/js arrays",
            xy=(0.50, 0.76),
            wh=(0.40, 0.10),
        ),
        Box(
            id="stackloop",
            text="Loop over local stacks (y)\nfor each stack:\nset start_cell index",
            xy=(0.50, 0.64),
            wh=(0.44, 0.12),
            fc="#F1F5F9",
        ),
        Box(
            id="interp",
            text="Build fine grid (z):\ninterpolate material params\nD, τ_sf, τ_φ, J_sd, β_c, …\n+ magnetization m̂(z)",
            xy=(0.26, 0.49),
            wh=(0.44, 0.14),
        ),
        Box(
            id="charge",
            text="Charge transient (ns, ne)\nImplicit backward-Euler\n+ laser source Q(z,t)\nCompute charge current Jc(z)",
            xy=(0.74, 0.49),
            wh=(0.44, 0.14),
        ),
        Box(
            id="split",
            text="Spin accumulation S(z,t)\n(Strang splitting)",
            xy=(0.50, 0.37),
            wh=(0.44, 0.09),
            fc="#FFF7ED",
        ),
        Box(
            id="diff1",
            text="Diffusion half-step\nImplicit solve for S\n+ interface coupling\n+ Neumann boundaries",
            xy=(0.22, 0.26),
            wh=(0.44, 0.12),
        ),
        Box(
            id="heun",
            text="Explicit Heun step\n∂S/∂t from:\n-∂J_drift/∂z, demag source\nprecession (S×m̂)\nspin-flip & dephasing",
            xy=(0.50, 0.20),
            wh=(0.44, 0.16),
        ),
        Box(
            id="diff2",
            text="Diffusion half-step\nImplicit solve for S",
            xy=(0.78, 0.26),
            wh=(0.44, 0.12),
        ),
        Box(
            id="coarse",
            text="Coarse average to microcells:\nsa_final(cell) = ⟨S⟩\nCompute Js = -β_c Jc m̂ - D ∂S/∂z\nFill ns/ne/jc/js outputs",
            xy=(0.50, 0.07),
            wh=(0.78, 0.14),
            fc="#F0FDFA",
        ),
        Box(
            id="torque",
            text="Convert sa_final → spin_torque field\n(used by field.cpp)\nspin_torque ∝ J_sd · sa_final",
            xy=(0.85, 0.07),
            wh=(0.28, 0.14),
            fc="#DCFCE7",
        ),
        Box(
            id="mpi",
            text="MPI gather (if MPICF)\nAllreduce per-cell arrays:\nspin_torque, sa_final,\nns/ne/jc/js",
            xy=(0.80, 0.76),
            wh=(0.34, 0.12),
            fc="#FFE4E6",
        ),
    ]

    edges = [
        Edge("init", "entry"),
        Edge("entry", "clear"),
        Edge("clear", "stackloop"),
        Edge("stackloop", "interp"),
        Edge("stackloop", "charge"),
        Edge("interp", "split"),
        Edge("charge", "split"),
        Edge("split", "diff1", label="dt/2"),
        Edge("diff1", "heun", label="dt"),
        Edge("heun", "diff2", label="dt/2"),
        Edge("diff2", "coarse"),
        Edge("coarse", "torque"),
        Edge("coarse", "mpi", label="after stacks"),
    ]

    return boxes, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--out",
        default="spincurrents-flowchart.png",
        help="Output PNG path",
    )
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    boxes, edges = build_flow()
    boxes_by_id = {b.id: b for b in boxes}

    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Draw edges first so arrows go under boxes.
    for e in edges:
        draw_edge(ax, boxes_by_id, e)

    for b in boxes:
        draw_box(ax, b)

    ax.text(
        0.01,
        0.99,
        "VAMPIRE spincurrents.cpp (Lepadatu-style 1D transient solver) — algorithm flow",
        ha="left",
        va="top",
        fontsize=13,
        weight="bold",
        color="#0B1220",
    )

    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

