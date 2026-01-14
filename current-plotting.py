#!/usr/bin/env python3
"""
1D spin currents Spin Transport Data Plotter

This script reads output from the VAMPIRE spincurrents.cpp 1D spin currents solver
and plots the z-profiles of charge and spin quantities.

Comparison with literature:
- 1D spin currents, S. PRB 2017: Transient spin accumulation in multilayers
- Remy et al.: Ultrafast demagnetization dynamics

Usage:
    python current-plotting.py [options]

The script looks for files in spin-acc/sc1d_data_* and creates
profile plots of ns, ne, Jc, and Js vs z.
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
import argparse

# Physical constants for unit conversion
MU_B = 9.274e-24  # Bohr magneton (J/T)
E_CHARGE = 1.602e-19  # Elementary charge (C)
HBAR = 1.055e-34  # Reduced Planck constant (J·s)

def read_sc1d_data(filename):
    """
    Read 1D spin currents output file.
    
    Columns: pos_x, pos_y, pos_z, ns, ne, Jc, Js_x, Js_y, Js_z
    
    Returns: numpy array of shape (N, 9) or empty array if no data
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 9:
                print(f"Warning: skipping malformed line in {filename}")
                continue
            try:
                data.append([float(x) for x in parts[:9]])
            except ValueError as e:
                print(f"Warning: could not parse line in {filename}: {e}")
                continue
    return np.array(data) if data else np.array([])


def extract_z_profiles(data):
    """
    Extract z-profiles by averaging over (x,y) positions at same z.
    
    Returns: dict with keys 'z', 'ns', 'ne', 'jc', 'js_x', 'js_y', 'js_z', 'js_mag'
    """
    if data.size == 0:
        return None
    
    # Get unique z values (with tolerance for floating point)
    z_vals = data[:, 2]
    z_unique = np.unique(np.round(z_vals, 6))
    
    # Average quantities at each z
    n_z = len(z_unique)
    result = {
        'z': np.zeros(n_z),
        'ns': np.zeros(n_z),
        'ne': np.zeros(n_z),
        'jc': np.zeros(n_z),
        'js_x': np.zeros(n_z),
        'js_y': np.zeros(n_z),
        'js_z': np.zeros(n_z),
        'js_mag': np.zeros(n_z)
    }
    
    for i, z in enumerate(z_unique):
        mask = np.abs(z_vals - z) < 1e-6
        result['z'][i] = z
        result['ns'][i] = np.mean(data[mask, 3])
        result['ne'][i] = np.mean(data[mask, 4])
        result['jc'][i] = np.mean(data[mask, 5])
        result['js_x'][i] = np.mean(data[mask, 6])
        result['js_y'][i] = np.mean(data[mask, 7])
        result['js_z'][i] = np.mean(data[mask, 8])
        result['js_mag'][i] = np.sqrt(result['js_x'][i]**2 + 
                                       result['js_y'][i]**2 + 
                                       result['js_z'][i]**2)
    
    # Sort by z
    sort_idx = np.argsort(result['z'])
    for key in result:
        result[key] = result[key][sort_idx]
    
    return result


def get_time_from_filename(fname):
    """Extract timestep number from filename."""
    match = re.search(r'sc1d_data_(\d+)', fname)
    return int(match.group(1)) if match else 0


def plot_profiles(files_to_plot, output_prefix='sc1d'):
    """
    Create profile plots for selected time steps.
    
    Following 1D spin currents PRB 2017 conventions:
    - ns: non-equilibrium spin density (related to ultrafast demagnetization)
    - ne: excess charge density
    - Jc: charge current density
    - Js: spin current density (vector)
    """
    if not files_to_plot:
        print("No 1D spin currents data files found in spin-acc/")
        return

    # Sort files numerically by timestep
    files_to_plot = sorted(files_to_plot, key=get_time_from_filename)
    
    # Select representative time steps (first, middle, last)
    if len(files_to_plot) > 5:
        indices = np.linspace(0, len(files_to_plot)-1, 5, dtype=int)
        selected_files = [files_to_plot[i] for i in indices]
    else:
        selected_files = files_to_plot
    
    print(f"Found {len(files_to_plot)} files, plotting {len(selected_files)}:")
    for f in selected_files:
        print(f"  - {os.path.basename(f)}")

    # Create figure with 2x2 subplots
    fig, axs = plt.subplots(2, 2, figsize=(14, 11))
    axs = axs.flatten()
    
    # Color map for different time steps
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(selected_files)))

    for idx, fname in enumerate(selected_files):
        data = read_sc1d_data(fname)
        profiles = extract_z_profiles(data)
        
        if profiles is None:
            print(f"Warning: no valid data in {fname}")
            continue
        
        # z is already in Angstroms from the output file
        z_nm = profiles['z'] * 0.1  # Å to nm
        
        timestep = get_time_from_filename(fname)
        label = f't = {timestep}'
        color = colors[idx]
        
        # Plot ns (non-equilibrium charge density)
        axs[0].plot(z_nm, profiles['ns'], label=label, color=color, linewidth=1.5)
        
        # Plot ne (excess charge density)  
        axs[1].plot(z_nm, profiles['ne'], label=label, color=color, linewidth=1.5)
        
        # Plot Jc (charge current density)
        axs[2].plot(z_nm, profiles['jc'], label=label, color=color, linewidth=1.5)
        
        # Plot |Js| (spin current magnitude)
        axs[3].plot(z_nm, profiles['js_mag'], label=label, color=color, linewidth=1.5)

    # Configure subplots
    titles = [
        r'Non-equilibrium Charge Density $n_s$',
        r'Excess Charge Density $n_e$',
        r'Charge Current Density $J_c$',
        r'Spin Current Magnitude $|J_s|$'
    ]
    ylabels = [
        r'$n_s$ (m$^{-3}$)',
        r'$n_e$ (m$^{-3}$)',
        r'$J_c$ (A/m$^2$)',
        r'$|J_s|$ (A/m$^2$)'
    ]
    
    for i, ax in enumerate(axs):
        ax.set_xlabel('z (nm)', fontsize=11)
        ax.set_ylabel(ylabels[i], fontsize=11)
        ax.set_title(titles[i], fontsize=12)
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(-2, 2))

    plt.tight_layout()
    
    # Save figure
    outfile = f'{output_prefix}_profiles.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"\nProfile plot saved to: {outfile}")
    
    plt.show()


def plot_spin_current_components(files_to_plot, output_prefix='sc1d'):
    """
    Plot all three components of spin current separately.
    
    This is useful for comparing with 1D spin currents's results which show
    spin polarization direction effects.
    """
    if not files_to_plot:
        print("No files found.")
        return
    
    # Use only the last file (final state)
    files_to_plot = sorted(files_to_plot, key=get_time_from_filename)
    fname = files_to_plot[-1]
    
    data = read_sc1d_data(fname)
    profiles = extract_z_profiles(data)
    
    if profiles is None:
        print("No valid data.")
        return
    
    z_nm = profiles['z'] * 0.1
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(z_nm, profiles['js_x'], 'r-', label=r'$J_s^x$', linewidth=2)
    ax.plot(z_nm, profiles['js_y'], 'g-', label=r'$J_s^y$', linewidth=2)
    ax.plot(z_nm, profiles['js_z'], 'b-', label=r'$J_s^z$', linewidth=2)
    ax.plot(z_nm, profiles['js_mag'], 'k--', label=r'$|J_s|$', linewidth=1.5, alpha=0.7)
    
    ax.set_xlabel('z (nm)', fontsize=12)
    ax.set_ylabel(r'$J_s$ (A/m$^2$)', fontsize=12)
    ax.set_title(f'Spin Current Components (t = {get_time_from_filename(fname)})', fontsize=13)
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(-2, 2))
    
    plt.tight_layout()
    outfile = f'{output_prefix}_js_components.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Spin current components plot saved to: {outfile}")
    plt.show()


def analyze_dynamics(files_to_plot):
    """
    Analyze time evolution of integrated quantities.
    
    This mimics the analysis in 1D spin currents papers where they track
    the total spin current at interfaces vs time.
    """
    if not files_to_plot:
        print("No files found.")
        return
    
    files_to_plot = sorted(files_to_plot, key=get_time_from_filename)
    
    timesteps = []
    total_js = []
    peak_ns = []
    
    for fname in files_to_plot:
        data = read_sc1d_data(fname)
        profiles = extract_z_profiles(data)
        
        if profiles is None:
            continue
        
        timesteps.append(get_time_from_filename(fname))
        
        # Integrate |Js| over z (simplified - assumes uniform dz)
        dz = np.diff(profiles['z'])
        if len(dz) > 0:
            js_integrated = np.sum(profiles['js_mag'][:-1] * dz)
        else:
            js_integrated = 0
        total_js.append(js_integrated)
        
        # Peak ns value
        peak_ns.append(np.max(np.abs(profiles['ns'])))
    
    if not timesteps:
        print("No valid data for dynamics analysis.")
        return
    
    timesteps = np.array(timesteps)
    total_js = np.array(total_js)
    peak_ns = np.array(peak_ns)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    ax1.plot(timesteps, total_js, 'b-o', markersize=4)
    ax1.set_ylabel(r'$\int |J_s| dz$ (A/m)', fontsize=11)
    ax1.set_title('Time Evolution of Spin Transport', fontsize=12)
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(timesteps, peak_ns, 'r-o', markersize=4)
    ax2.set_xlabel('Time step', fontsize=11)
    ax2.set_ylabel(r'Peak $|n_s|$ (m$^{-3}$)', fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sc1d_dynamics.png', dpi=150, bbox_inches='tight')
    print("Dynamics plot saved to: sc1d_dynamics.png")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot 1D spin currents spin transport data')
    parser.add_argument('--dir', default='spin-acc', help='Directory containing data files')
    parser.add_argument('--components', action='store_true', help='Plot Js components separately')
    parser.add_argument('--dynamics', action='store_true', help='Analyze time evolution')
    args = parser.parse_args()
    
    # Find data files
    pattern = os.path.join(args.dir, 'sc1d_data_*')
    files = glob.glob(pattern)
    
    if not files:
        print(f"No files matching {pattern}")
        print("Make sure the simulation has run with 1D spin currents solver enabled.")
    else:
        print(f"Found {len(files)} 1D spin currents data files\n")
        
        # Main profile plots
        plot_profiles(files)
        
        # Optional: Js components
        if args.components:
            plot_spin_current_components(files)
        
        # Optional: dynamics analysis
        if args.dynamics:
            analyze_dynamics(files)
