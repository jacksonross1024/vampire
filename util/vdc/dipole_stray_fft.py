#!/usr/bin/env python3
"""
Periodic (tiled flake) magnetostatic stray field above a multilayer, via 2D FFT.

Input A — **cells-*.txt** (same as cells-pub-ori.py):
  - Skip header line, tab-separated.
  - Use rows with cellz == 0 (single microcell layer in-plane).
  - Per physical layer L in {0..3}: five columns starting at index 3+5*L:
      mx, my, mz, |m|, n_atoms
  - Magnetic moment per cell (Bohr magneton number) along i in {x,y,z}:
        mu_i = m_i * |m| * n_atoms
  - In-plane positions: x_nm = 0.1 * cellx, y_nm = 0.1 * celly (same as cells-pub-ori.py).
  - Cell surface area for normalisation: normalise (nm^2), default 16.0 as in that script.

Input B — **avg-magnetisation-spatial-grid.txt** (from avg-cells.py):
  - Comment lines ``#`` then a header row: x, y, L1_mx..L1_mz, L2_*, L3_*, L4_* (14 columns).
  - Each Lk_m* is ``2.98 * (component moment) / normalise`` with normalise=16 in avg-cells.py,
    where component moment is ``m_i * |m| * n_atoms`` in μ_B number (same as cells).
  - Recover μ (μ_B number): ``mu = Lk_* * normalise / 2.98`` (overridable via --avg-scale and
    the same --normalise as cells surface area).
  - x, y are already in nm. For the m_x,m_y,m_z PNG, plotted values are the file columns
    (same units as the grid file), not raw unit-cell m vectors.

Layer heights (nm): z0 = 0, z1 = z0 + dz, z2 = z1 + dz, z3 = z2 + dz with dz = 6.54.

Stray field **B** (magnetic flux density in **tesla**) is a single 2D map in one vacuum plane:
z_obs = z_top + height_above_nm, where z_top is the top magnetisation layer height (layer
index 3). Values **Hx, Hy, Hz** in the output NPZ are **B_x, B_y, B_z** in vacuum
(**B** = μ₀ **H**, with **H** the magnetostatic field strength from the k-space kernel below).
There is no separate “dipole layer index”: only this height sets where **B** is evaluated.
Optional PNGs: (1) **m_x, m_y, m_z** from the input file (unitless; ``--plot-m-layer``);
(2) **B_x, B_y, B_z** in **tesla** at ``z_obs`` from the FFT (use ``--no-plot-b`` to skip).

The k-space kernel is the usual 2D periodic magnetostatic continuation (Laplace above the
stack): for each Fourier mode k = (kx, ky), k = |k|, decay = exp(-k * |z_obs - z_l|), and
surface dipole densities px, py, pz (A), the vacuum field is

  Hx = -mu0 * sum_l decay_l * ( (kx^2*px + kx*ky*py)/k^2 + (kx*pz)/k )
  Hy = -mu0 * sum_l decay_l * ( (kx*ky*px + ky^2*py)/k^2 + (ky*pz)/k )
  Hz = -mu0 * sum_l decay_l * ( (kx*px + ky*py)/k + k*pz )

with the k=0 mode set to zero (no net monopole / stray field DC component).

**Coarse-graining (optional):** ``--coarse-dx-nm`` / ``--coarse-dy-nm`` bin microcells
before building the FFT grid. For each bin, **net moment** is the **sum** of per-cell
``mu_i = m_i*|m|*n_atoms`` (μ_B). **Surface dipole density** uses that total moment divided
by the **sum of microcell areas** in the bin (``n_cells_in_bin * --normalise`` nm²), i.e.
areal average of μ_B/nm², consistent with the FFT expecting p = (total A·m²) / (physical area).

Requires: numpy, matplotlib
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
except ImportError as e:  # pragma: no cover
    plt = None  # type: ignore
    LinearSegmentedColormap = None  # type: ignore
    Normalize = None  # type: ignore
    _MPL_IMPORT_ERROR = e
else:
    _MPL_IMPORT_ERROR = None

# Physical constants (SI)
MU0 = 4.0 * np.pi * 1e-7  # H/m
MU_B = 9.2740100783e-24  # J/T = A*m^2

# avg-cells.py: Lk_* = AVG_MOMENT_SCALE * (m_i * |m| * n_atoms) / AVG_GRID_NORMALISE
DEFAULT_AVG_MOMENT_SCALE = 2.98
DEFAULT_AVG_GRID_NORMALISE = 16.0


def _peek_first_non_comment_line(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            return s
    return ""


def detect_input_format(path: str) -> str:
    """
    'cells': raw cells-*.txt (≥23 tab columns on first data row).
    'avg_grid': avg-magnetisation-spatial-grid.txt (14 cols; header x, L1_mx, …).
    """
    s = _peek_first_non_comment_line(path)
    if not s:
        raise ValueError(f"{path}: empty or only comments")
    parts = s.split("\t")
    if len(parts) >= 23:
        return "cells"
    if len(parts) == 14:
        if parts[0].strip().lower() == "x" and "L1_mx" in s:
            return "avg_grid"
        try:
            float(parts[0])
            float(parts[1])
            return "avg_grid"
        except ValueError:
            pass
    return "cells"


def load_avg_magnetisation_spatial_grid(
    path: str,
    x_width: float | None = None,
    y_width: float | None = None,
    num_layers: int = 4,
    *,
    avg_moment_scale: float = DEFAULT_AVG_MOMENT_SCALE,
    avg_grid_normalise: float = DEFAULT_AVG_GRID_NORMALISE,
):
    """
    Load avg-cells.py output: x, y (nm), L1_mx..L4_mz (scaled as in avg-cells.py).

    Recovers μ in Bohr magneton number: mu_i = file_value * avg_grid_normalise / avg_moment_scale.

    Returns:
      x_nm, y_nm, mu, m_raw — same shapes as load_cells_layers.
      m_raw holds file columns Lk_mx, Lk_my, Lk_mz (for plotting).
    """
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    i0 = 0
    while i0 < len(lines) and (
        not lines[i0].strip() or lines[i0].lstrip().startswith("#")
    ):
        i0 += 1
    if i0 >= len(lines):
        raise ValueError(f"{path}: no header")
    header_line = lines[i0].strip()
    if "L1_mx" not in header_line:
        raise ValueError(f"{path}: expected avg grid header with L1_mx (line {i0 + 1})")
    data_start = i0 + 1
    arr = np.loadtxt(path, dtype=np.float64, delimiter="\t", skiprows=data_start)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] != 14:
        raise ValueError(f"{path}: expected 14 columns, got {arr.shape[1]}")

    xw = float("inf") if x_width is None else float(x_width)
    yw = float("inf") if y_width is None else float(y_width)

    inv_scale = float(avg_grid_normalise) / float(avg_moment_scale)

    xs_nm: list[float] = []
    ys_nm: list[float] = []
    mus: list[list[np.ndarray]] = []
    m_raws: list[list[np.ndarray]] = []

    for row in range(arr.shape[0]):
        x_nm = float(arr[row, 0])
        y_nm = float(arr[row, 1])
        if x_nm > xw or y_nm > yw:
            continue

        layer_mu: list[np.ndarray] = []
        layer_m_raw: list[np.ndarray] = []
        for L in range(num_layers):
            o = 2 + 3 * L
            fx, fy, fz = arr[row, o], arr[row, o + 1], arr[row, o + 2]
            layer_m_raw.append(np.array([fx, fy, fz], dtype=np.float64))
            layer_mu.append(
                np.array([fx * inv_scale, fy * inv_scale, fz * inv_scale], dtype=np.float64)
            )

        xs_nm.append(x_nm)
        ys_nm.append(y_nm)
        mus.append(layer_mu)
        m_raws.append(layer_m_raw)

    if not xs_nm:
        raise RuntimeError(f"No rows left in {path} after x/y bounds filter.")

    x_nm = np.asarray(xs_nm, dtype=np.float64)
    y_nm = np.asarray(ys_nm, dtype=np.float64)
    mu = np.asarray(mus, dtype=np.float64)
    m_raw = np.asarray(m_raws, dtype=np.float64)
    return x_nm, y_nm, mu, m_raw


def load_cells_layers(
    path: str,
    x_width: float | None = None,
    y_width: float | None = None,
    num_layers: int = 4,
):
    """
    Load cell positions and per-layer moments (mu_B number) from cells-*.txt.

    Returns:
      x_nm, y_nm: 1D arrays of positions for used cells
      mu: (n_cells, num_layers, 3) moment components in units of Bohr magneton number
      m_raw: (n_cells, num_layers, 3) raw mx, my, mz from the file (before * |m| * n_atoms)
    """
    mag_data = np.genfromtxt(path, dtype=float, unpack=True, delimiter="\t", skip_header=1)
    ncols = mag_data.shape[0]
    need = 3 + 5 * num_layers
    if ncols < need:
        raise ValueError(f"{path}: expected at least {need} columns, got {ncols}")

    xw = float("inf") if x_width is None else float(x_width)
    yw = float("inf") if y_width is None else float(y_width)

    xs_nm: list[float] = []
    ys_nm: list[float] = []
    mus: list[list[np.ndarray]] = []
    m_raws: list[list[np.ndarray]] = []

    count = 0
    nrows = mag_data.shape[1]
    for i in range(nrows):
        if mag_data[2][i] != 0.0:
            count += 1
            continue
        x_nm = 0.1 * mag_data[0][i]
        y_nm = 0.1 * mag_data[1][i]
        if x_nm > xw or y_nm > yw:
            count += 1
            continue

        layer_mu: list[np.ndarray] = []
        layer_m_raw: list[np.ndarray] = []
        for L in range(num_layers):
            b = 3 + 5 * L
            mx, my, mz = mag_data[b][i], mag_data[b + 1][i], mag_data[b + 2][i]
            ml = mag_data[b + 3][i]
            nat = mag_data[b + 4][i]
            mu_vec = np.array(
                [mx * ml * nat, my * ml * nat, mz * ml * nat],
                dtype=np.float64,
            )
            layer_mu.append(mu_vec)
            layer_m_raw.append(np.array([mx, my, mz], dtype=np.float64))

        xs_nm.append(x_nm)
        ys_nm.append(y_nm)
        mus.append(layer_mu)
        m_raws.append(layer_m_raw)
        count += 1

    if not xs_nm:
        raise RuntimeError(f"No valid rows in {path} (cellz==0 and within bounds).")

    x_nm = np.asarray(xs_nm, dtype=np.float64)
    y_nm = np.asarray(ys_nm, dtype=np.float64)
    mu = np.asarray(mus, dtype=np.float64)  # (n, num_layers, 3)
    m_raw = np.asarray(m_raws, dtype=np.float64)  # (n, num_layers, 3)
    return x_nm, y_nm, mu, m_raw


def coarse_grain_cell_bins(
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    mu: np.ndarray,
    m_raw: np.ndarray,
    coarse_dx_nm: float,
    coarse_dy_nm: float,
    cell_area_nm2: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    int,
    int,
    float,
    float,
    np.ndarray,
    np.ndarray,
]:
    """
    Bin microcells into coarse rectangles; sum μ (μ_B) per bin; average raw m for plotting.

    Dipole density for each bin uses total moment / (n_cells * cell_area_nm2).

    Returns
    -------
    x_nm, y_nm : (nx*ny,) bin centres
    mu : (nx*ny, num_layers, 3) summed μ_B per bin
    m_raw : (nx*ny, num_layers, 3) mean raw m in bin (for PNG)
    area_nm2 : (nx*ny,) total microcell area per bin (0 if empty)
    nx, ny, dx_nm, dy_nm
    ix, iy : (nx*ny,) integer indices for placing into (nx, ny) grids
    """
    x_nm = np.asarray(x_nm, dtype=np.float64)
    y_nm = np.asarray(y_nm, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    m_raw = np.asarray(m_raw, dtype=np.float64)
    n, num_layers, three = mu.shape
    if three != 3 or m_raw.shape != (n, num_layers, 3):
        raise ValueError("mu / m_raw shape mismatch")
    if coarse_dx_nm <= 0 or coarse_dy_nm <= 0:
        raise ValueError("coarse bin sizes must be positive")
    x0 = float(np.min(x_nm))
    y0 = float(np.min(y_nm))
    ix = np.floor((x_nm - x0) / float(coarse_dx_nm)).astype(np.int64, copy=False)
    iy = np.floor((y_nm - y0) / float(coarse_dy_nm)).astype(np.int64, copy=False)
    pairs = np.stack([ix, iy], axis=1)
    uniq, inv, counts = np.unique(
        pairs, axis=0, return_inverse=True, return_counts=True
    )
    n_bins = uniq.shape[0]
    nx = int(uniq[:, 0].max() + 1)
    ny = int(uniq[:, 1].max() + 1)

    mu_sum = np.zeros((n_bins, num_layers, 3), dtype=np.float64)
    m_acc = np.zeros((n_bins, num_layers, 3), dtype=np.float64)
    for L in range(num_layers):
        for c in range(3):
            mu_sum[:, L, c] = np.bincount(
                inv, weights=mu[:, L, c], minlength=n_bins
            )
            m_acc[:, L, c] = np.bincount(
                inv, weights=m_raw[:, L, c], minlength=n_bins
            )
    with np.errstate(invalid="ignore"):
        m_mean = m_acc / counts[:, np.newaxis, np.newaxis]

    mu_dense = np.zeros((nx, ny, num_layers, 3), dtype=np.float64)
    m_raw_dense = np.zeros((nx, ny, num_layers, 3), dtype=np.float64)
    count_dense = np.zeros((nx, ny), dtype=np.int64)
    for b in range(n_bins):
        bi, bj = int(uniq[b, 0]), int(uniq[b, 1])
        mu_dense[bi, bj] = mu_sum[b]
        m_raw_dense[bi, bj] = m_mean[b]
        count_dense[bi, bj] = int(counts[b])

    xc = x0 + (np.arange(nx, dtype=np.float64) + 0.5) * float(coarse_dx_nm)
    yc = y0 + (np.arange(ny, dtype=np.float64) + 0.5) * float(coarse_dy_nm)
    II, JJ = np.meshgrid(
        np.arange(nx, dtype=np.int64),
        np.arange(ny, dtype=np.int64),
        indexing="ij",
    )
    ix_flat = II.ravel()
    iy_flat = JJ.ravel()
    x_out = xc[ix_flat]
    y_out = yc[iy_flat]
    mu_out = mu_dense.reshape(-1, num_layers, 3)
    m_out = m_raw_dense.reshape(-1, num_layers, 3)
    area_nm2 = (count_dense.ravel().astype(np.float64)) * float(cell_area_nm2)
    return (
        x_out,
        y_out,
        mu_out,
        m_out,
        area_nm2,
        nx,
        ny,
        float(coarse_dx_nm),
        float(coarse_dy_nm),
        ix_flat,
        iy_flat,
    )


def infer_grid(
    x_nm: np.ndarray,
    y_nm: np.ndarray,
):
    """Infer regular grid from unique coordinates; returns strides (nm), nx, ny, ix, iy."""
    ux = np.unique(np.round(x_nm, 9))
    uy = np.unique(np.round(y_nm, 9))
    if len(ux) < 2 or len(uy) < 2:
        raise RuntimeError("Need at least two distinct x and y coordinates to infer grid.")
    dx_nm = float(np.min(np.diff(np.sort(ux))))
    dy_nm = float(np.min(np.diff(np.sort(uy))))
    ix = np.rint(x_nm / dx_nm).astype(np.int64)
    iy = np.rint(y_nm / dy_nm).astype(np.int64)
    ix -= ix.min()
    iy -= iy.min()
    nx = int(ix.max() + 1)
    ny = int(iy.max() + 1)
    return dx_nm, dy_nm, nx, ny, ix, iy


def moments_to_surface_A(
    mu_layers: np.ndarray,
    area_nm2: float | np.ndarray,
):
    """
    mu_layers: (n, 3) or (num_layers, 3) moment in mu_B number per row.
    area_nm2: scalar, or (n,) total area in nm^2 per row (sum of microcell areas when
        coarse-grained). Rows with zero area get zero surface density.
    Returns px, py, pz each shape (n,) or (num_layers,) in A (dipole moment per unit area).
    """
    mu_layers = np.asarray(mu_layers, dtype=np.float64)
    mu_SI = mu_layers * MU_B
    a = np.asarray(area_nm2, dtype=np.float64)
    if a.ndim == 0 or a.size == 1:
        area_m2 = float(a.ravel()[0]) * 1e-18
        return mu_SI / area_m2
    if a.shape[0] != mu_SI.shape[0]:
        raise ValueError(
            f"area_nm2 length {a.shape[0]} != mu_layers rows {mu_SI.shape[0]}"
        )
    area_m2 = a * 1e-18
    out = np.zeros_like(mu_SI)
    mask = area_m2 > 0.0
    out[mask] = mu_SI[mask] / area_m2[mask, np.newaxis]
    return out


def stray_field_fft_layers(
    px_grids: list[np.ndarray],
    py_grids: list[np.ndarray],
    pz_grids: list[np.ndarray],
    dx_m: float,
    dy_m: float,
    z_layers_m: np.ndarray,
    z_obs_m: float,
):
    """
    Sum stray-field contribution from each layer (same nx, ny).

    Returns Hx, Hy, Hz from the magnetostatic kernel in **A/m**; caller converts to
    **B** in **tesla** with B = μ₀ H for NPZ output.
    """
    nx, ny = px_grids[0].shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx_m)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy_m)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    k2 = KX * KX + KY * KY
    k = np.sqrt(k2)
    k_safe = np.maximum(k, 1e-30)

    Hx_k = np.zeros((nx, ny), dtype=np.complex128)
    Hy_k = np.zeros((nx, ny), dtype=np.complex128)
    Hz_k = np.zeros((nx, ny), dtype=np.complex128)

    for px, py, pz, zl in zip(px_grids, py_grids, pz_grids, z_layers_m):
        px_k = np.fft.fft2(px)
        py_k = np.fft.fft2(py)
        pz_k = np.fft.fft2(pz)
        dz = float(z_obs_m - zl)
        decay = np.exp(-k_safe * abs(dz))

        inv_k = 1.0 / k_safe
        inv_k2 = 1.0 / np.maximum(k2, 1e-60)

        # Avoid NaN at k=0: contributions are multiplied by k in numerator for some terms
        Hx_k += decay * (
            -MU0 * (KX * KX * px_k + KX * KY * py_k) * inv_k2
            - MU0 * (KX * pz_k) * inv_k
        )
        Hy_k += decay * (
            -MU0 * (KX * KY * px_k + KY * KY * py_k) * inv_k2
            - MU0 * (KY * pz_k) * inv_k
        )
        Hz_k += decay * (
            -MU0 * (KX * px_k + KY * py_k) * inv_k
            - MU0 * k_safe * pz_k
        )

    # Enforce no DC stray-field component
    Hx_k[0, 0] = 0.0
    Hy_k[0, 0] = 0.0
    Hz_k[0, 0] = 0.0

    Hx = np.real(np.fft.ifft2(Hx_k))
    Hy = np.real(np.fft.ifft2(Hy_k))
    Hz = np.real(np.fft.ifft2(Hz_k))
    return Hx, Hy, Hz


def _magnetization_cmap_norm(vmin: float, vmax: float):
    """Colormap style aligned with cells-pub-ori.py (blue–white–purple)."""
    assert LinearSegmentedColormap is not None and Normalize is not None
    mag_min = vmin
    mag_max = vmax
    if mag_min < 0.0:
        rng = mag_max - mag_min
        if abs(mag_min) / rng > 0.1:
            nodes = [
                0.0,
                abs(mag_min) / rng,
                abs(mag_min) / rng + 0.333 * mag_max / rng,
                abs(mag_min) / rng + 0.667 * mag_max / rng,
                1.0,
            ]
            cmap = LinearSegmentedColormap.from_list(
                "mycolors",
                list(zip(nodes, ["#81B1CB", "#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"])),
            )
            norm = Normalize(vmin=mag_min, vmax=mag_max)
        else:
            nodes = [0.0, 0.333, 0.667, 1.0]
            cmap = LinearSegmentedColormap.from_list(
                "mycolors",
                list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"])),
            )
            norm = Normalize(vmin=0.0, vmax=mag_max)
    else:
        nodes = [0.0, 0.333, 0.667, 1.0]
        cmap = LinearSegmentedColormap.from_list(
            "mycolors",
            list(zip(nodes, ["#f5f0f7", "#B5AFCF", "#7A5FA9", "#614199"])),
        )
        norm = Normalize(vmin=mag_min, vmax=mag_max)
    return cmap, norm


def plot_m_components_png(
    mx_grid: np.ndarray,
    my_grid: np.ndarray,
    mz_grid: np.ndarray,
    x_max_nm: float,
    y_max_nm: float,
    out_path: str,
    title_suffix: str = "",
):
    """1x3 figure: m_x, m_y, m_z (same layout as cells-pub-ori layer panels)."""
    if plt is None:
        raise RuntimeError(f"matplotlib is required for plotting: {_MPL_IMPORT_ERROR}")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    extent = (0.0, x_max_nm, 0.0, y_max_nm)
    labels = (r"$m_x$", r"$m_y$", r"$m_z$")
    grids = (mx_grid, my_grid, mz_grid)
    for ax, lab, grid in zip(axes, labels, grids):
        vmin = float(np.nanmin(grid))
        vmax = float(np.nanmax(grid))
        if abs(vmax - vmin) < 1e-30:
            vmax = vmin + 1.0
        cmap, norm = _magnetization_cmap_norm(vmin, vmax)
        im = ax.imshow(
            grid.T,
            extent=extent,
            origin="lower",
            cmap=cmap,
            norm=norm,
            interpolation="gaussian",
            aspect="equal",
        )
        ax.set_xlabel("x (nm)")
        ax.set_ylabel("y (nm)")
        ax.set_title(lab + (f" {title_suffix}" if title_suffix else ""))
        cb = fig.colorbar(im, ax=ax, shrink=0.8, aspect=12)
        cb.set_label("(unitless)")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def plot_b_components_png(
    bx_grid: np.ndarray,
    by_grid: np.ndarray,
    bz_grid: np.ndarray,
    x_max_nm: float,
    y_max_nm: float,
    out_path: str,
    title_suffix: str = "",
):
    """1x3 figure: B_x, B_y, B_z (tesla) at the observation plane."""
    if plt is None:
        raise RuntimeError(f"matplotlib is required for plotting: {_MPL_IMPORT_ERROR}")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    extent = (0.0, x_max_nm, 0.0, y_max_nm)
    labels = (r"$B_x$", r"$B_y$", r"$B_z$")
    grids = (bx_grid, by_grid, bz_grid)
    for ax, lab, grid in zip(axes, labels, grids):
        zf = np.asarray(grid, dtype=np.float64)[np.isfinite(grid)]
        if zf.size:
            vmax = float(np.nanpercentile(np.abs(zf), 99.0))
        else:
            vmax = 1.0
        if vmax <= 0 or not np.isfinite(vmax):
            vmax = 1e-12
        im = ax.imshow(
            grid.T,
            extent=extent,
            origin="lower",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            interpolation="gaussian",
            aspect="equal",
        )
        ax.set_xlabel("x (nm)")
        ax.set_ylabel("y (nm)")
        ax.set_title(lab + (f" {title_suffix}" if title_suffix else ""))
        cb = fig.colorbar(im, ax=ax, shrink=0.8, aspect=12)
        cb.set_label("T")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="FFT stray field above periodic cells multilayer.")
    p.add_argument(
        "cells",
        nargs="?",
        default="cells-00000000.txt",
        help=(
            "Input: cells-*.txt or avg-magnetisation-spatial-grid.txt (see --input-format)"
        ),
    )
    p.add_argument(
        "--input-format",
        choices=("auto", "cells", "avg_grid"),
        default="auto",
        help="auto: detect from columns (14 + L1_mx vs ≥23); cells: raw VAMPIRE cells file",
    )
    p.add_argument(
        "--avg-moment-scale",
        type=float,
        default=DEFAULT_AVG_MOMENT_SCALE,
        metavar="S",
        help=(
            "avg_grid input only: same S as avg-cells.py (Lk_* = S * mu_muB / normalise); "
            f"default {DEFAULT_AVG_MOMENT_SCALE}"
        ),
    )
    p.add_argument(
        "--height-above-top-nm",
        type=float,
        default=10.0,
        help=(
            "Vacuum observation height: distance (nm) above the top magnetisation layer "
            "(index 3) where Bx,By,Bz are evaluated (output in teslas) — one horizontal plane only "
            "(default: 10)"
        ),
    )
    p.add_argument(
        "--layer-spacing-nm",
        type=float,
        default=6.54,
        help="Vertical spacing between adjacent layers in nm (default: 6.54)",
    )
    p.add_argument(
        "--normalise",
        type=float,
        default=16.0,
        help="Cell surface area in nm^2 (default: 16, same as cells-pub-ori.py)",
    )
    p.add_argument(
        "--x-width",
        type=float,
        default=856.0,
        help="Max x position in nm (same as cells-pub-ori.py default)",
    )
    p.add_argument(
        "--y-width",
        type=float,
        default=1188.0,
        help="Max y position in nm (same as cells-pub-ori.py default)",
    )
    p.add_argument(
        "--out-prefix",
        type=str,
        default="stray_fft",
        help="Prefix for output .npz (default: stray_fft)",
    )
    p.add_argument(
        "--plot-m-layer",
        type=int,
        default=0,
        help=(
            "For the m_x,m_y,m_z figure only: which magnetisation layer (0..3) in the cells "
            "file to plot. Does not set the dipole-field height; use --height-above-top-nm for that."
        ),
    )
    p.add_argument(
        "--no-plot-m",
        action="store_true",
        help="Skip writing the 1x3 m_x,m_y,m_z PNG",
    )
    p.add_argument(
        "--no-plot-b",
        action="store_true",
        help="Skip writing the 1x3 B_x,B_y,B_z PNG (stray flux density in teslas at z_obs)",
    )
    p.add_argument(
        "--coarse-dx-nm",
        type=float,
        default=None,
        metavar="DX",
        help=(
            "Bin microcells in x before FFT (nm). Total μ_B per bin is summed; "
            "dipole density uses sum(μ)/(n_cells_in_bin * --normalise) nm². "
            "Default --coarse-dy-nm: same as DX"
        ),
    )
    p.add_argument(
        "--coarse-dy-nm",
        type=float,
        default=None,
        metavar="DY",
        help="Bin microcells in y (nm); default: same as --coarse-dx-nm",
    )
    args = p.parse_args(argv)

    path = args.cells
    if not os.path.isfile(path):
        print(f"Error: file not found: {path}", file=sys.stderr)
        return 1

    fmt = args.input_format
    if fmt == "auto":
        fmt = detect_input_format(path)
    if fmt == "avg_grid":
        x_nm, y_nm, mu, m_raw = load_avg_magnetisation_spatial_grid(
            path,
            x_width=args.x_width,
            y_width=args.y_width,
            num_layers=4,
            avg_moment_scale=args.avg_moment_scale,
            avg_grid_normalise=float(args.normalise),
        )
    else:
        x_nm, y_nm, mu, m_raw = load_cells_layers(
            path,
            x_width=args.x_width,
            y_width=args.y_width,
            num_layers=4,
        )

    use_coarse = (
        args.coarse_dx_nm is not None or args.coarse_dy_nm is not None
    )
    area_nm2_for_p: float | np.ndarray = float(args.normalise)
    if use_coarse:
        cdx = args.coarse_dx_nm
        cdy = args.coarse_dy_nm
        if cdx is None:
            cdx = cdy
        if cdy is None:
            cdy = cdx
        (
            x_nm,
            y_nm,
            mu,
            m_raw,
            area_nm2_for_p,
            nx,
            ny,
            dx_nm,
            dy_nm,
            ix,
            iy,
        ) = coarse_grain_cell_bins(
            x_nm,
            y_nm,
            mu,
            m_raw,
            float(cdx),
            float(cdy),
            float(args.normalise),
        )
    else:
        dx_nm, dy_nm, nx, ny, ix, iy = infer_grid(x_nm, y_nm)
    dx_m = dx_nm * 1e-9
    dy_m = dy_nm * 1e-9

    num_layers = mu.shape[1]
    z_layers_nm = np.arange(num_layers, dtype=np.float64) * args.layer_spacing_nm
    z_top_nm = z_layers_nm[-1]
    z_obs_nm = z_top_nm + args.height_above_top_nm
    z_layers_m = z_layers_nm * 1e-9
    z_obs_m = z_obs_nm * 1e-9

    px_grids: list[np.ndarray] = []
    py_grids: list[np.ndarray] = []
    pz_grids: list[np.ndarray] = []

    for L in range(num_layers):
        px = np.zeros((nx, ny), dtype=np.float64)
        py = np.zeros((nx, ny), dtype=np.float64)
        pz = np.zeros((nx, ny), dtype=np.float64)
        p_A = moments_to_surface_A(mu[:, L, :], area_nm2_for_p)
        px[ix, iy] = p_A[:, 0]
        py[ix, iy] = p_A[:, 1]
        pz[ix, iy] = p_A[:, 2]
        px_grids.append(px)
        py_grids.append(py)
        pz_grids.append(pz)

    Hx, Hy, Hz = stray_field_fft_layers(
        px_grids,
        py_grids,
        pz_grids,
        dx_m,
        dy_m,
        z_layers_m,
        z_obs_m,
    )
    # Kernel returns H (A/m); store B = μ₀ H in vacuum (tesla) under keys Hx, Hy, Hz.
    Hx = Hx * MU0
    Hy = Hy * MU0
    Hz = Hz * MU0

    stem = os.path.basename(path).replace(".txt", "")
    out_npz = f"{args.out_prefix}_H_{stem}.npz"
    np.savez(
        out_npz,
        Hx=Hx,
        Hy=Hy,
        Hz=Hz,
        x_nm=np.arange(nx, dtype=np.float64) * dx_nm,
        y_nm=np.arange(ny, dtype=np.float64) * dy_nm,
        z_obs_nm=z_obs_nm,
        z_layers_nm=z_layers_nm,
        dx_nm=dx_nm,
        dy_nm=dy_nm,
        normalise_nm2=args.normalise,
        coarse_dx_nm=np.array(float("nan") if not use_coarse else float(dx_nm)),
        coarse_dy_nm=np.array(float("nan") if not use_coarse else float(dy_nm)),
        cells_file=os.path.abspath(path),
        input_format=np.array(fmt),
        avg_moment_scale=np.array(float(args.avg_moment_scale)),
        Hx_Hy_Hz_units_Tesla=np.array(True),
    )
    print(f"Wrote {out_npz}")
    print(
        f"  grid {nx}x{ny}, dx={dx_nm:.6g} nm, dy={dy_nm:.6g} nm, "
        f"z_layers_nm={z_layers_nm}, z_obs_nm={z_obs_nm:.6g}"
    )
    if use_coarse:
        print(
            "  coarse-grain: summed μ_B per bin; p uses total_area = n_cells * normalise "
            f"(nm²); bin size {dx_nm:g} × {dy_nm:g} nm"
        )
    print(
        f"  B range (T): x [{Hx.min():.4g}, {Hx.max():.4g}], "
        f"y [{Hy.min():.4g}, {Hy.max():.4g}], z [{Hz.min():.4g}, {Hz.max():.4g}]"
    )

    x_max_nm = float(nx) * dx_nm
    y_max_nm = float(ny) * dy_nm
    if not args.no_plot_b:
        b_png = (
            f"{args.out_prefix}_B_components_{os.path.basename(path).replace('.txt', '')}.png"
        )
        try:
            plot_b_components_png(
                Hx,
                Hy,
                Hz,
                x_max_nm,
                y_max_nm,
                b_png,
                title_suffix=f"(z_obs={z_obs_nm:g} nm)",
            )
            print(f"Wrote {b_png}")
        except RuntimeError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    if not args.no_plot_m:
        if args.plot_m_layer < 0 or args.plot_m_layer >= num_layers:
            print(
                f"Error: --plot-m-layer must be in 0..{num_layers - 1}",
                file=sys.stderr,
            )
            return 1
        mx_g = np.zeros((nx, ny), dtype=np.float64)
        my_g = np.zeros((nx, ny), dtype=np.float64)
        mz_g = np.zeros((nx, ny), dtype=np.float64)
        Lp = args.plot_m_layer
        mx_g[ix, iy] = m_raw[:, Lp, 0]
        my_g[ix, iy] = m_raw[:, Lp, 1]
        mz_g[ix, iy] = m_raw[:, Lp, 2]
        png_path = f"{args.out_prefix}_m_components_{os.path.basename(path).replace('.txt', '')}.png"
        try:
            plot_m_components_png(
                mx_g,
                my_g,
                mz_g,
                x_max_nm,
                y_max_nm,
                png_path,
                title_suffix=f"(layer {Lp})",
            )
            print(f"Wrote {png_path}")
        except RuntimeError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
