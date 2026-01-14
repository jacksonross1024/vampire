#!/usr/bin/env python3
"""
Generate a simplified, publication-style flowchart for the 1D transient
spin-transport model (Lepadatu Eq. 1–5 style) with an added demag-driven
spin-accumulation source term.

Output PNG(s) intended for quick communication of the physics, not code detail.

No external deps beyond matplotlib.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Optional, Tuple, List

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


@dataclass
class Node:
    key: str
    title: str
    body: str
    xy: Tuple[float, float]
    wh: Tuple[float, float]
    fc: str
    ec: str = "#0F172A"


@dataclass
class Link:
    a: str
    b: str
    label: Optional[str] = None


def _box(ax, n: Node):
    x, y = n.xy
    w, h = n.wh
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.6,
        edgecolor=n.ec,
        facecolor=n.fc,
    )
    ax.add_patch(patch)

    ax.text(
        x,
        y + h * 0.20,
        n.title,
        ha="center",
        va="center",
        fontsize=12.5,
        weight="bold",
        color="#0B1220",
        family="DejaVu Sans",
    )
    ax.text(
        x,
        y - h * 0.12,
        n.body,
        ha="center",
        va="center",
        fontsize=10.3,
        color="#0B1220",
        family="DejaVu Sans",
        linespacing=1.25,
    )


def _anchor(n: Node, which: str) -> Tuple[float, float]:
    x, y = n.xy
    w, h = n.wh
    if which == "top":
        return (x, y + h / 2)
    if which == "bottom":
        return (x, y - h / 2)
    if which == "left":
        return (x - w / 2, y)
    if which == "right":
        return (x + w / 2, y)
    raise ValueError(which)


def _arrow(ax, na: Node, nb: Node, label: Optional[str]):
    # Always draw vertical arrows top→bottom in this diagram.
    p1 = _anchor(na, "bottom")
    p2 = _anchor(nb, "top")
    arr = FancyArrowPatch(
        p1,
        p2,
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=1.35,
        color="#0F172A",
        connectionstyle="arc3,rad=0.0",
    )
    ax.add_patch(arr)
    if label:
        ax.text(
            (p1[0] + p2[0]) / 2,
            (p1[1] + p2[1]) / 2,
            label,
            ha="center",
            va="center",
            fontsize=9.5,
            color="#334155",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9),
        )


def build_nodes() -> Tuple[List[Node], List[Link]]:
    # This is intentionally “physics first”, aligned to a Lepadatu-style
    # drift–diffusion formulation (often presented as Eq. 1–5).
    nodes = [
        Node(
            key="inputs",
            title="Inputs / Materials",
            body=(
                "Layer stack + microcells (1D along z)\n"
                "D(z), β(z), τsf(z), τφ(z), Jsd(z)\n"
                "m̂(z,t) from atomistic LLG\n"
                "Laser / electrical excitation"
            ),
            xy=(0.5, 0.88),
            wh=(0.78, 0.16),
            fc="#EEF2FF",
        ),
        Node(
            key="charge",
            title="Charge sector (Eq. 1–3)",
            body=(
                "Solve transient charge densities\n"
                "n_s(z,t), n_e(z,t) (implicit)\n"
                "→ charge current J_c(z,t)\n"
                "including source Q(z,t) (laser)"
            ),
            xy=(0.5, 0.66),
            wh=(0.78, 0.16),
            fc="#ECFEFF",
        ),
        Node(
            key="spin_pde",
            title="Spin accumulation PDE (Eq. 4)",
            body=(
                "Evolve spin accumulation S(z,t)\n"
                "∂S/∂t = ∂/∂z(D∂S/∂z)  − ∂J_drift/∂z\n"
                "        − S/τsf − S⊥/τφ  − ω(S×m̂)\n"
                "+ demag-driven source:  −χ · (d|m|/dt) · m̂_ref"
            ),
            xy=(0.5, 0.44),
            wh=(0.78, 0.20),
            fc="#FFF7ED",
        ),
        Node(
            key="numerics",
            title="Time stepping",
            body=(
                "Strang splitting:\n"
                "Diffusion (implicit, dt/2)\n"
                "→ Explicit Heun (dt)\n"
                "→ Diffusion (implicit, dt/2)\n"
                "Interface coupling + Neumann outer BCs"
            ),
            xy=(0.5, 0.24),
            wh=(0.78, 0.16),
            fc="#F1F5F9",
        ),
        Node(
            key="outputs",
            title="Outputs (Eq. 5 + torque map)",
            body=(
                "Spin current:  J_s = −β J_c m̂ − D ∂S/∂z\n"
                "Store: n_s, n_e, J_c, J_s, S\n"
                "Torque field to atoms:\n"
                "H_ST ∝ Jsd · S   (added to LLG field)"
            ),
            xy=(0.5, 0.06),
            wh=(0.78, 0.16),
            fc="#DCFCE7",
        ),
    ]

    links = [
        Link("inputs", "charge"),
        Link("charge", "spin_pde", label="Jc drives drift term"),
        Link("spin_pde", "numerics"),
        Link("numerics", "outputs"),
    ]
    return nodes, links


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="spincurrents-flowchart-lite.png")
    ap.add_argument("--dpi", type=int, default=260)
    args = ap.parse_args()

    nodes, links = build_nodes()
    by_key = {n.key: n for n in nodes}

    fig = plt.figure(figsize=(14.5, 7.8))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Title
    ax.text(
        0.02,
        0.985,
        "1D transient spin transport (Lepadatu-style Eq. 1–5) + demag-driven source",
        ha="left",
        va="top",
        fontsize=16,
        weight="bold",
        color="#0B1220",
        family="DejaVu Sans",
    )
    ax.text(
        0.02,
        0.955,
        "Simplified physics flow (matches spincurrents.cpp high-level structure)",
        ha="left",
        va="top",
        fontsize=11.5,
        color="#334155",
        family="DejaVu Sans",
    )

    # Links underneath
    for lk in links:
        _arrow(ax, by_key[lk.a], by_key[lk.b], lk.label)

    # Nodes on top
    for n in nodes:
        _box(ax, n)

    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

