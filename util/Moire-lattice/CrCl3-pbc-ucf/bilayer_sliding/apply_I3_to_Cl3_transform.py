#!/usr/bin/env python3
"""
I3 → Cl3 bilayer_sliding maps in *normalised lattice* space.

Crystal lattices (lookup uses physical Å = u · a):
    a_Cl3 = 5.94 Å,  a_I3 = 6.93 Å

Inter isotropic fit (shared reduced u = r_xy / a only; Δz / c never enters):
    Sample Kim Same/Diff at r = u · a_Cl3
    Sample I3        at r = R(θ)·(u)·a_I3 + (dx, dy)
    Fit only on |u| ≤ RCUT/a_Cl3 = 5.2/5.94  (Cl3 inter DFT horizon)
        I3 ≈ α · Same + β · Diff + c
    k_amp = std(αS+βD) / std(I3)  on that mask

Cl3 tables are written on the shared reader grid (±A0,±A1) but *indexed as
Cl3 physical Å* so match_inter(dx,dy) with Cl3 bonds is correct:
    u = (x,y) / a_Cl3
    Inter J  : Kim (αS+βD+c) for |u|≤u_cut;
               k_amp · I3(u) for |u|>u_cut  (I3 long-range, distance-scaled)
    Inter DMI: k_amp · I3(u) everywhere (incl. long-range; no |r| wipe)
    Intra    : k_amp · I3; sliding axes × (a_Cl3/a_I3)  (Cl3 moiré slide units)

I3 long-range source: bilayer_sliding_fullrange/ (interpolate_inter_v2 --r-cut 0).
Optional rot/dx/dy from kim_fig6 sweep (fine alignment only; a fixed to crystal).
"""

from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator

N = 100
A0, A1 = 7.276, 7.402  # reader half-widths (Å); Cl3 values live at Cl physical (x,y)
A_CL3 = 5.94
A_I3 = 6.93  # crystal; NOT the old fit free-a (~7.28)
RCUT = 5.2  # Cl3 inter DFT horizon (Å); fit mask |u|≤RCUT/a_Cl3

ROOT = Path(__file__).resolve().parents[1]
I3DIR = ROOT.parent / "CrI3-pbc-ucf" / "bilayer_sliding"
I3FULL = ROOT.parent / "CrI3-pbc-ucf" / "bilayer_sliding_fullrange"
CLDIR = Path(__file__).resolve().parent
META = CLDIR / "kim_fig6"


def load_r(p: Path) -> np.ndarray:
    return np.loadtxt(p).reshape(N, N)


def save_r(p: Path, arr: np.ndarray) -> None:
    with open(p, "w") as f:
        for v in np.ravel(arr):
            f.write(f"{v:.4f}\n")


def reflect_J(J, cr):
    if cr == 1:
        return J
    if cr == 2:
        return np.flipud(J)
    if cr == 3:
        return np.fliplr(J)
    if cr == 4:
        return np.flipud(np.fliplr(J))
    raise ValueError(cr)


def reflect_D(Dx, Dy, Dz, cr):
    """Species mirrors matching interpolate_inter_v2 (axial-like under y/x flip)."""
    if cr == 1:
        return Dx, Dy, Dz
    if cr == 2:
        return -np.flipud(Dx), np.flipud(Dy), -np.flipud(Dz)
    if cr == 3:
        return -np.fliplr(Dx), np.fliplr(Dy), -np.fliplr(Dz)
    if cr == 4:
        return (
            np.flipud(np.fliplr(Dx)),
            np.flipud(np.fliplr(Dy)),
            np.flipud(np.fliplr(Dz)),
        )
    raise ValueError(cr)


def add_offset_preserve_zero(L, c, soft=None):
    """Apply additive c without pushing near-zero L across zero."""
    L = np.asarray(L, dtype=float)
    if soft is None:
        soft = max(abs(float(c)) * 0.5, 1e-4)
    w = np.abs(L) / (np.abs(L) + soft)
    J = L + float(c) * w
    flipped = (L * J < 0) & (np.abs(L) < soft)
    return np.where(flipped, 0.0, J)


def soft_approach_zero(D, floor=1e-4):
    D = np.asarray(D, dtype=float)
    return np.where(np.abs(D) < floor, 0.0, D)


def load_spatial_params():
    """rot/dx/dy from sweep; a_I3 is always crystal A_I3 (ignore sweep best_a)."""
    npz = META / "I3_a_rot_xy_sweep.npz"
    if npz.is_file():
        z = np.load(npz)
        return float(z["best_ang"]), float(z["best_dx"]), float(z["best_dy"])
    return -2.0, 0.30, -0.10


def make_rgi(field: np.ndarray, x_coords: np.ndarray, y_coords: np.ndarray):
    """RegularGridInterpolator on ascending (y, x); field[i,j] with i=0 at high y."""
    return RegularGridInterpolator(
        (y_coords, x_coords),
        field[::-1, :],
        bounds_error=False,
        fill_value=np.nan,
    )


def ensure_i3_fullrange() -> Path:
    """Build uncapped I3 inter tables if missing (needed for long-range Cl3)."""
    need = [
        "Cr1_inter_map_2.txt",
        "Cr1_Dx_inter_map_2_avg.txt",
        "Cr1_Dy_inter_map_2_avg.txt",
        "Cr1_Dz_inter_map_2_avg.txt",
    ]
    if I3FULL.is_dir() and all((I3FULL / f).is_file() for f in need):
        return I3FULL
    import subprocess
    import sys

    script = (
        ROOT.parent / "CrI3" / "bilayer_sliding_2" / "interpolate_inter_v2.py"
    )
    I3FULL.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(script),
        "--r-cut",
        "0",
        "--out",
        str(I3FULL),
        "--src",
        str(ROOT.parent / "CrI3" / "bilayer_sliding"),
    ]
    print("building I3 full-range inter maps:", " ".join(cmd))
    subprocess.check_call(cmd)
    return I3FULL


def main():
    META.mkdir(exist_ok=True)
    ROT_DEG, DX, DY = load_spatial_params()
    SCALE_LAT = A_CL3 / A_I3
    u_cut_Cl = RCUT / A_CL3
    i3_src = ensure_i3_fullrange()

    jj, ii = np.meshgrid(np.arange(N), np.arange(N))
    # Reader grid physical Å — for Cl3 tables this *is* Cl3 bond (dx,dy)
    x_A = jj * 0.02 * A0 - A0
    y_A = (N - 1 - ii) * 0.02 * A1 - A1
    r_A = np.sqrt(x_A**2 + y_A**2)
    x_coords = x_A[0, :]
    y_coords = y_A[::-1, 0]

    ng = 201
    us = np.linspace(-1.0, 1.0, ng)
    vs = np.linspace(1.0, -1.0, ng)  # row0 = +u_y
    UU, VV = np.meshgrid(us, vs)
    RR_u = np.sqrt(UU**2 + VV**2)

    th = np.deg2rad(ROT_DEG)
    cth, sth = np.cos(th), np.sin(th)
    Ur = cth * UU - sth * VV
    Vr = sth * UU + cth * VV

    print(
        f"lattice: a_I3={A_I3}  a_Cl3={A_CL3}  SCALE_LAT={SCALE_LAT:.6f}"
    )
    print(
        f"align: rot={ROT_DEG:.2f}°  dx={DX:.3f} Å  dy={DY:.3f} Å  "
        f"(fine alignment only)"
    )
    print(
        f"fit disk |u|≤{u_cut_Cl:.4f} (= {RCUT}/{A_CL3}); "
        f"I3 at same u → r≤{u_cut_Cl * A_I3:.3f} Å"
    )
    print(f"I3 full-range maps: {i3_src}")

    S0 = load_r(CLDIR / "Kim_same_inter_map.txt")
    D0 = load_r(CLDIR / "Kim_diff_inter_map.txt")
    I0 = load_r(i3_src / "Cr1_inter_map_2.txt")
    interp_S = make_rgi(S0, x_coords, y_coords)
    interp_D = make_rgi(D0, x_coords, y_coords)
    interp_I = make_rgi(I0, x_coords, y_coords)

    # Shared reduced-u sampling (standardised to each lattice)
    pts_Cl = np.column_stack([(VV * A_CL3).ravel(), (UU * A_CL3).ravel()])
    S_u = interp_S(pts_Cl).reshape(UU.shape)
    D_u = interp_D(pts_Cl).reshape(UU.shape)

    pts_I = np.column_stack(
        [(Vr * A_I3 + DY).ravel(), (Ur * A_I3 + DX).ravel()]
    )
    I_u = interp_I(pts_I).reshape(UU.shape)

    fit_mask = (
        np.isfinite(S_u)
        & np.isfinite(D_u)
        & np.isfinite(I_u)
        & (RR_u <= u_cut_Cl)
        & (RR_u <= 1.0)
    )
    print(f"fit points: {fit_mask.sum()}")

    Sv, Dv, Iv = S_u[fit_mask], D_u[fit_mask], I_u[fit_mask]
    A = np.column_stack([Sv, Dv, np.ones_like(Sv)])
    coef, _, _, _ = np.linalg.lstsq(A, Iv, rcond=None)
    alpha, beta, c0 = [float(x) for x in coef]
    pred = alpha * Sv + beta * Dv + c0
    rmse = float(np.sqrt(np.mean((pred - Iv) ** 2)))
    corr = float(np.corrcoef(pred, Iv)[0, 1])
    print(f"I3 = {alpha:.6f}*Same + {beta:.6f}*Diff + {c0:.6f}")
    print(f"RMSE={rmse:.6f}  corr={corr:.6f}")

    L_u = alpha * S_u + beta * D_u
    I_from_Cl_u = add_offset_preserve_zero(L_u, c0)
    std_L = float(np.std(L_u[fit_mask]))
    std_I = float(np.std(I_u[fit_mask]))
    k_amp = std_L / std_I if std_I > 1e-12 else 1.0
    print(f"k_amp = std(αS+βD)/std(I3) = {k_amp:.6f}")

    def sample_i3_on_cl_grid(field_i3: np.ndarray, fill=0.0) -> np.ndarray:
        """I3 field at reduced u of each Cl reader point (x_A,y_A)/a_Cl3."""
        interp = make_rgi(field_i3, x_coords, y_coords)
        u = x_A / A_CL3
        v = y_A / A_CL3
        ur = cth * u - sth * v
        vr = sth * u + cth * v
        pts = np.column_stack(
            [(vr * A_I3 + DY).ravel(), (ur * A_I3 + DX).ravel()]
        )
        out = interp(pts).reshape(N, N)
        return np.where(np.isfinite(out), out, fill)

    def kim_on_cl_grid(field_kim: np.ndarray, fill=0.0) -> np.ndarray:
        """Kim maps already live on reader Å = Cl physical lookup coords."""
        interp = make_rgi(field_kim, x_coords, y_coords)
        pts = np.column_stack([y_A.ravel(), x_A.ravel()])
        out = interp(pts).reshape(N, N)
        return np.where(np.isfinite(out), out, fill)

    # --- Inter J on Cl physical grid ---
    S_cl = kim_on_cl_grid(S0)
    D_cl = kim_on_cl_grid(D0)
    L_cl = alpha * S_cl + beta * D_cl
    J_kim = add_offset_preserve_zero(L_cl, c0)
    J_i3 = sample_i3_on_cl_grid(I0)
    u_r = r_A / A_CL3
    J1_cl = soft_approach_zero(
        np.where(u_r <= u_cut_Cl, J_kim, k_amp * J_i3)
    )
    n_long = int(np.sum(u_r > u_cut_Cl))
    n_long_nz = int(np.sum((u_r > u_cut_Cl) & (np.abs(J1_cl) > 1e-8)))
    print(
        f"inter J: Kim inside |u|≤{u_cut_Cl:.4f}; "
        f"k_amp·I3 outside ({n_long} cells, {n_long_nz} nonzero)"
    )

    # --- Inter DMI: k_amp · I3(u) everywhere (long-range kept) ---
    Dx1 = load_r(i3_src / "Cr1_Dx_inter_map_2_avg.txt")
    Dy1 = load_r(i3_src / "Cr1_Dy_inter_map_2_avg.txt")
    Dz1 = load_r(i3_src / "Cr1_Dz_inter_map_2_avg.txt")
    Dx1_cl = soft_approach_zero(k_amp * sample_i3_on_cl_grid(Dx1))
    Dy1_cl = soft_approach_zero(k_amp * sample_i3_on_cl_grid(Dy1))
    Dz1_cl = soft_approach_zero(k_amp * sample_i3_on_cl_grid(Dz1))
    print(
        f"inter DMI: k_amp·I3 on Cl grid; "
        f"|D|>0 cells={(np.abs(Dx1_cl) + np.abs(Dy1_cl) + np.abs(Dz1_cl) > 1e-8).sum()}"
    )

    for cr in [1, 2, 3, 4]:
        save_r(CLDIR / f"Cr{cr}_inter_map_2.txt", reflect_J(J1_cl, cr))
        Dx, Dy, Dz = reflect_D(Dx1_cl, Dy1_cl, Dz1_cl, cr)
        save_r(CLDIR / f"Cr{cr}_Dx_inter_map_2_avg.txt", Dx)
        save_r(CLDIR / f"Cr{cr}_Dy_inter_map_2_avg.txt", Dy)
        save_r(CLDIR / f"Cr{cr}_Dz_inter_map_2_avg.txt", Dz)
        print(f"wrote Cr{cr} inter (Cl3 physical Å, normalised-u mapped)")

    # Intra: sliding axes in Cl3 lattice units; values × k_amp
    def transform_intra(src: Path, dst: Path):
        d = np.loadtxt(src)
        assert d.shape[0] == 12 * N * N
        out = d.copy()
        out[:, 0] *= SCALE_LAT
        out[:, 1] *= SCALE_LAT
        out[:, 2] = soft_approach_zero(out[:, 2] * k_amp)
        out[:, 3] = soft_approach_zero(out[:, 3] * k_amp)
        out[:, 4] = soft_approach_zero(out[:, 4] * k_amp)
        out[:, 5] = soft_approach_zero(out[:, 5] * k_amp)
        np.savetxt(dst, out, fmt="%.4f", delimiter="\t")

    for cr in [1, 2, 3, 4]:
        transform_intra(
            I3DIR / f"Cr{cr}_intra_data.txt",
            CLDIR / f"Cr{cr}_intra_data.txt",
        )
        print(f"wrote Cr{cr}_intra_data.txt (axes×{SCALE_LAT:.4f}, ×k_amp)")

    from symmetrize_sliding_maps import symmetrize_inter, symmetrize_intra, verify

    print("symmetrizing sliding maps (C3 + species)...")
    symmetrize_inter()
    symmetrize_intra()
    verify()

    np.savez(
        META / "fit_Cl_to_I3_constant_amp.npz",
        alpha=alpha,
        beta=beta,
        c=c0,
        k_amp=k_amp,
        rmse=rmse,
        corr=corr,
        a_I3=A_I3,
        a_Cl3=A_CL3,
        rot_deg=ROT_DEG,
        dx=DX,
        dy=DY,
        rcut=RCUT,
        u_cut=u_cut_Cl,
        scale_lat=SCALE_LAT,
    )
    with open(META / "fit_Cl_to_I3_constant_amp.txt", "w") as f:
        f.write("# Normalised-u fit: I3(u) = alpha*Same(u) + beta*Diff(u) + c\n")
        f.write(
            f"# a_I3={A_I3} a_Cl3={A_CL3} rot={ROT_DEG} deg "
            f"dx={DX} dy={DY}; fit |u|<={u_cut_Cl} (RCUT={RCUT}Å / a_Cl3)\n"
        )
        f.write(f"alpha {alpha:.10g}\n")
        f.write(f"beta {beta:.10g}\n")
        f.write(f"c {c0:.10g}\n")
        f.write(f"k_amp {k_amp:.10g}\n")
        f.write(f"RMSE {rmse:.10g}\n")
        f.write(f"corr {corr:.10g}\n")
        f.write(f"a_I3 {A_I3:.10g}\n")
        f.write(f"a_Cl3 {A_CL3:.10g}\n")
        f.write(f"rot_deg {ROT_DEG:.10g}\n")
        f.write(f"dx {DX:.10g}\n")
        f.write(f"dy {DY:.10g}\n")
        f.write(f"u_cut {u_cut_Cl:.10g}\n")
        f.write(
            "# Inter J: Kim αS+βD+c for |u|≤u_cut; k_amp·I3(u) beyond "
            "(Cl3 physical table)\n"
        )
        f.write(
            "# Inter DMI + Intra: k_amp·I3 in normalised space; "
            f"intra axes × {SCALE_LAT:.10g}\n"
        )

    # ----- plots -----
    import matplotlib.pyplot as plt

    def plot_2x4(cl_maps, i3_maps, labels, title, outfile):
        fig, axs = plt.subplots(2, 4, figsize=(14, 7))
        for j, lab in enumerate(labels):
            cl, i3m = cl_maps[j], i3_maps[j]
            vmax = max(
                np.percentile(np.abs(cl), 99),
                np.percentile(np.abs(i3m), 99),
                1e-3,
            )
            im0 = axs[0, j].imshow(
                cl, origin="upper", cmap="RdBu_r", vmin=-vmax, vmax=vmax
            )
            axs[0, j].set_title(f"Cl3 {lab}", fontsize=11)
            axs[0, j].set_xticks([])
            axs[0, j].set_yticks([])
            plt.colorbar(im0, ax=axs[0, j], fraction=0.046, pad=0.02)
            im1 = axs[1, j].imshow(
                i3m, origin="upper", cmap="RdBu_r", vmin=-vmax, vmax=vmax
            )
            axs[1, j].set_title(f"I3 {lab}", fontsize=11)
            axs[1, j].set_xticks([])
            axs[1, j].set_yticks([])
            plt.colorbar(im1, ax=axs[1, j], fraction=0.046, pad=0.02)
        axs[0, 0].set_ylabel("Cl3", fontsize=12)
        axs[1, 0].set_ylabel("I3", fontsize=12)
        fig.suptitle(title, fontsize=12)
        fig.tight_layout()
        fig.savefig(outfile, dpi=160, bbox_inches="tight")
        print("wrote", outfile)

    labels = ["J (isotropic)", "Dx", "Dy", "Dz"]
    inter_files = [
        "Cr{}_inter_map_2.txt",
        "Cr{}_Dx_inter_map_2_avg.txt",
        "Cr{}_Dy_inter_map_2_avg.txt",
        "Cr{}_Dz_inter_map_2_avg.txt",
    ]

    for cr in [1, 2, 3, 4]:
        files = [f.format(cr) for f in inter_files]
        # Compare Cl3 (Cl Å) to I3 fullrange sampled via same reflect on I3 grid
        plot_2x4(
            [load_r(CLDIR / f) for f in files],
            [load_r(i3_src / f) for f in files],
            labels,
            f"INTER Cr{cr}: Cl3 (Cl-physical, u-mapped) vs I3 full-range\n"
            f"I3 = {alpha:.3f}·Same + {beta:.3f}·Diff + {c0:.3f}  "
            f"fit |u|≤{u_cut_Cl:.3f}; long-range J/DMI = k_amp·I3(u)  "
            f"k={k_amp:.3f}  a_I3={A_I3}, a_Cl3={A_CL3}",
            META / f"Cl3_vs_I3_INTER_Cr{cr}_J_DMI_2x4.png",
        )

    def intra_block(path, block=0):
        d = np.loadtxt(path)
        sl = slice(block * N * N, (block + 1) * N * N)
        b = d[sl]
        return [b[:, k].reshape(N, N) for k in (2, 3, 4, 5)]

    intra_shells = (
        ("1NN", 0),
        ("2NN", 3),
        ("3NN", 9),
    )
    for cr in [1, 2, 3, 4]:
        i3_path = I3DIR / f"Cr{cr}_intra_data.txt"
        cl_path = CLDIR / f"Cr{cr}_intra_data.txt"
        for shell, block in intra_shells:
            bi = np.loadtxt(i3_path)[block * N * N, :2]
            bc = np.loadtxt(cl_path)[block * N * N, :2]
            thb = float(np.degrees(np.arctan2(bc[1], bc[0])))
            plot_2x4(
                intra_block(cl_path, block),
                intra_block(i3_path, block),
                labels,
                f"INTRA Cr{cr} {shell} (block {block}, θ≈{thb:.0f}°): "
                f"Cl3 (×k={k_amp:.3f}, axes×{SCALE_LAT:.3f}) vs I3\n"
                f"bond I3=({bi[0]:.3f},{bi[1]:.3f})  Cl3=({bc[0]:.3f},{bc[1]:.3f})",
                META / f"Cl3_vs_I3_INTRA_Cr{cr}_{shell}_J_DMI_2x4.png",
            )
        legacy = META / f"Cl3_vs_I3_INTRA_Cr{cr}_J_DMI_2x4.png"
        src1 = META / f"Cl3_vs_I3_INTRA_Cr{cr}_1NN_J_DMI_2x4.png"
        if src1.is_file():
            import shutil

            shutil.copy2(src1, legacy)

    fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.0))
    vmax = max(np.nanpercentile(np.abs(I_u[fit_mask]), 99), 0.05)
    extent = [-1.0, 1.0, -1.0, 1.0]

    def show_u(ax, data, title, xlabel):
        disp = np.where(fit_mask, data, np.nan)
        im = ax.imshow(
            disp,
            origin="upper",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            extent=extent,
            aspect="equal",
        )
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("y/a")
        circ = plt.Circle((0, 0), u_cut_Cl, fill=False, ec="k", ls="--", lw=0.8)
        ax.add_patch(circ)
        plt.colorbar(im, ax=ax, fraction=0.046)
        return im

    show_u(
        axs[0],
        I_u,
        f"I3(u)  a={A_I3}, rot={ROT_DEG:.1f}°",
        "x/a (shared u)",
    )
    show_u(
        axs[1],
        I_from_Cl_u,
        f"αS+βD+c at a_Cl3={A_CL3}",
        "x/a_Cl3",
    )
    show_u(axs[2], I_u - I_from_Cl_u, "residual (I3 − fit)", "x/a")
    fig.suptitle(
        f"Normalised-u fit |u|≤{u_cut_Cl:.3f}: α={alpha:.3f}, β={beta:.3f}, "
        f"c={c0:.3f}, RMSE={rmse:.3f}, k_amp={k_amp:.3f}\n"
        f"dashed = Cl3 inter DFT horizon (r={RCUT}Å / a_Cl3)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(META / "Cl_to_I3_constant_fit.png", dpi=150, bbox_inches="tight")
    print("wrote fit diagnostic")


if __name__ == "__main__":
    main()
