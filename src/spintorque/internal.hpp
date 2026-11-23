#ifndef SPINTORQUE_INTERNAL_H_
#define SPINTORQUE_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

#include <string>
#include <vector>

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// spin-torque implementation. These functions should not be accessed
// outside of the spin-torque code.
//---------------------------------------------------------------------
namespace st{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Shared variables used for the spin torque calculation
      //-----------------------------------------------------------------------------
      extern bool enabled; // enable spin torque calculation
      extern bool TMRenable;

      extern bool fbc;
      //extern double micro_cell_size; /// lateral size of spin torque microcells
      extern std::vector<double> micro_cell_size;
      extern std::string microcell_decomp_type;
      extern double micro_cell_thickness; /// thickness of spin torque microcells (atomistic)

      extern int num_local_atoms; /// number of local atoms (ignores halo atoms in parallel simulation)
      extern int remove_nm;
      extern int current_direction; /// direction for current x->0, y->1, z->2
      //   std::vector< std::vector< micro_cell_t > > stack;
      extern std::vector<int> atom_st_index; // mc which atom belongs to
      extern std::vector<int> cell_stack_index; //stack each cell belongs to
      extern std::vector<double> x_field_array; // arrays to store atomic spin torque field
      extern std::vector<double> y_field_array;
      extern std::vector<double> z_field_array;

      extern int num_stacks_y;  // total number of stacks
      extern int num_stacks_x;  // total number of stacks
      extern int num_x_stacks; // number of stacks in x
      extern int num_y_stacks; // number of stack in y
      extern int num_microcells_per_stack; // number of microcells per stack

      extern int config_file_counter; // spin torque config file counter

      extern int stx;//=0; // indices for x,y,z in the spin torque coordinate system (default z)
      extern int sty;//=1;
      extern int stz;//=2;
      extern int free_layer;       /// index of free layer in magnetic tunnel junction
      extern int reference_layer;  /// index of reference layer in magnetic tunnel junction

      extern double je; // current (C/s)
      extern double initial_beta;
      extern double initial_theta;
      extern double rel_angle;
      extern int ST_output_rate;
      extern std::string output_torque_data;
      // spin-currents-1d-output. coarse: unit-cell bins. atomic: field windows. verbose: every fine cell.
      extern int sc1d_output_mode;
      const int sc1d_out_coarse = 0;
      const int sc1d_out_atomic = 1;
      const int sc1d_out_verbose = 2;

      extern std::vector<double> initial_m;
      extern std::vector<double> stack_init_mag;
      extern std::vector<double> init_stack_mag;

      extern std::vector<int> stack_index_y; // start of stack in microcell arrays
      extern std::vector<int> stack_index_x; // start of stack in microcell arrays
      extern std::vector<std::vector<int> > cell_index_x;

      //STT
      extern std::vector<double> beta_cond; /// spin polarisation (conductivity)
      extern std::vector<double> beta_diff; /// spin polarisation (diffusion)
      extern std::vector<double> sa_infinity; /// intrinsic spin accumulation
      extern std::vector<double> lambda_sdl; /// spin diffusion length
      extern std::vector<double> diffusion; /// spin diffusion length
      extern std::vector<double> sd_exchange; /// spin diffusion length
      extern std::vector<double> lambda_phi; /// transverse spin dephasing length (m)
      extern std::vector<double> chi_demag;  /// longitudinal demag source on S, (C/m^3) per muB
      extern std::vector<double> seebeck_coefficient; /// Seebeck coefficient S (V/K) for charge Seebeck effect
      extern std::vector<double> beta_c_seebeck; /// open-circuit Seebeck spin factor, set at init
      extern std::vector<double> superdiffusion_polarization; /// hot-electron polarization, set at init
      extern std::vector<double> conductivity; /// electrical conductivity σ (S/m); <0 Einstein from D, Te; 0 insulator
      extern std::vector<double> a; /// spin diffusion length
      extern std::vector<double> b; /// spin diffusion length

      //SOT
      extern std::vector<double> sot_beta_cond; /// spin polarisation (conductivity)
      extern std::vector<double> sot_beta_diff; /// spin polarisation (diffusion)
      extern std::vector<double> sot_sa_infinity; /// intrinsic spin accumulation
      extern std::vector<double> sot_lambda_sdl; /// spin diffusion length
      extern std::vector<double> sot_diffusion; /// spin diffusion length
      extern std::vector<double> sot_sd_exchange; /// spin diffusion length
      extern std::vector<double> sot_a; /// spin diffusion length
      extern std::vector<double> sot_b; /// spin diffusion length
      extern std::vector<double> spin_acc_sign;
      extern bool sot_sa;
      extern std::vector<bool> sot_sa_source;
      extern bool sot_check;

      // 1D spin accumulation solver (with optional demag-driven source)
      extern bool sc1d_enable;
      // Solve and write spin accumulation as usual, but do not add the
      // spin-torque field to the atom. Magnetisation still sources S.
      extern bool sc1d_calibrate;
      extern double sc1d_fine_dz;      // fine grid spacing in Angstroms (along st::internal::stz)
      // Prolongation between occupied fine cells: 0 step, 1 linear, 2 smoothstep.
      extern int sc1d_prolongation;
      const int sc1d_prolong_step = 0;
      const int sc1d_prolong_linear = 1;
      const int sc1d_prolong_smooth = 2;
      // Implicit diffusion half-step. 0 backward Euler, 1 BDF2.
      extern int sc1d_diffusion;
      const int sc1d_diff_be = 0;
      const int sc1d_diff_bdf2 = 1;
      extern bool sc1d_bdf2_primed;
      // Angstroms. Negative means the unit-cell height, resolved when the grid is built.
      // Step ignores it. Linear and smoothstep blend only inside this distance of a gap midpoint.
      extern double sc1d_prolongation_decay;
      // analytic broadens the material-pair sheet with the prolongation kernel.
      // diffuse averages those same constants by atomic population, then steps.
      extern int sc1d_interface;
      const int sc1d_interface_analytic = 0;
      const int sc1d_interface_diffuse = 1;
      extern int sc1d_spin_stride;     // update stride for spin accumulation (in LLG steps)
      extern int sc1d_charge_stride;   // update stride for charge/spin current fields (in LLG steps)
      extern double sc1d_temperature;  // user-defined temperature (K) placeholder for future thermal coupling
      extern unsigned long sc1d_relax_steps; // spin-current steps before the LLG is allowed to move spins

// Laser-driven charge transient parameters (Eqs. 1–4)
extern bool sc1d_laser_enable;
extern double sc1d_laser_Q0;                   // absorbed laser power density at z=0 (W/m^3), before time envelope
extern double sc1d_optical_absorption_length;  // optical absorption length d (m)
extern double sc1d_laser_wavelength;           // laser wavelength lambda (m)
extern double sc1d_laser_eta;                  // effective electrons excited per photon (dimensionless)
extern double sc1d_laser_t0;                   // pulse center time (s), relative to solver step counter
extern double sc1d_laser_fwhm;                 // pulse width (FWHM, s)
extern double sc1d_tau_s;                      // non-equilibrium electron lifetime tau_s (s)
extern double sc1d_tau_demag;                  // legacy input; target-channel lifetime, unused

// Superdiffusive transport parameters
extern bool sc1d_superdiffusive_enable;        // enable enhanced superdiffusive transport model
extern double sc1d_tau_e;                      // hot electron energy relaxation time (s) ~50-200 fs
extern double sc1d_v_e0;                       // base hot electron velocity (m/s), if 0 uses De/dopt
// spin-currents-1d-superdiffusive-type. 0 gradient, 1 constant speed, 2 lifetime decay.
// Gradient is ve = -De (dTe/dz) / Te times the post-peak speed factor.
// Constant and decay do not use De, so a mixing sheet cannot change them.
extern int sc1d_superdiffusive_type;
extern double sc1d_superdiffusive_speed;       // constant mode, input speed in m/s; the solver uses minus its absolute value
extern double sc1d_superdiffusive_lifetime;    // decay mode, s. |ve| goes from 1e6 to 1e5 over this time

// Temperature coupling parameters
extern bool sc1d_use_ltmp_temperatures;        // legacy input; temperature source is ltmp or the global TTM
extern double sc1d_reference_temperature;      // reference temperature for scaling (K)
extern bool sc1d_thermal_effects;              // transport constants may depend on Te, Tp (laws not filled in)

      // Seebeck effect parameters
      extern bool sc1d_seebeck_enable;               // enable charge Seebeck effect

      // Applied electrode potentials. Off keeps the series current at 0.
      extern bool sc1d_bias_enable;
      extern double sc1d_voltage_bottom;             // V at the first fine cell
      extern double sc1d_voltage_top;                // V at the last fine cell

      // Thermal gradient parameters
      extern bool sc1d_thermal_gradients_enable;      // enable thermal gradient calculation
      extern bool sc1d_thermal_gradients_initialised; // initialization flag
      extern bool sc1d_thermal_gradients_use_for_llg_fields; // use thermal gradients for LLG thermal fields (any program)

      // Thermal gradient arrays: [stack_idx * num_microcells_per_stack + cell_idx]
      extern std::vector<double> sc1d_Te;            // Electron temperature (K)
      extern std::vector<double> sc1d_Tp;            // Phonon temperature (K)
      extern std::vector<double> sc1d_sqrt_Te;       // sqrt(Te) for stability
      extern std::vector<double> sc1d_sqrt_Tp;       // sqrt(Tp) for stability

      // Material properties: [cell_idx] (same for all stacks, materials vary by z)
      extern std::vector<double> sc1d_Ce;            // Electron heat capacity (J/m³/K)
      extern std::vector<double> sc1d_Cp;            // Phonon heat capacity (J/m³/K)
      extern std::vector<double> sc1d_kappa_e;       // Electron thermal conductivity (J/s/m/K)
      extern std::vector<double> sc1d_kappa_p;       // Phonon thermal conductivity (J/s/m/K)
      extern std::vector<double> sc1d_G;             // Electron-phonon coupling (J/s/m³/K)
      extern std::vector<double> sc1d_T_Debye;       // Debye temperature (K)

      // Debye lookup table (shared, size ~24000)
      extern std::vector<double> sc1d_debye_table;

      // Atom mapping
      extern std::vector<int> sc1d_atom_cell_idx;    // Maps atom -> coarse cell index
      extern std::vector<bool> sc1d_atom_use_phonon; // true if atom couples to phonon temp

      // Fine-grid charge state per local stack (same nf as S)
      extern std::vector<double> sc1d_ns_fine;       // ns (C/m^3), size: local_stacks*nf
      extern std::vector<double> sc1d_ne_fine;       // ne (C/m^3), size: local_stacks*nf
      extern std::vector<double> sc1d_V_fine;        // electrostatic potential (V), size: local_stacks*nf
      extern std::vector<double> sc1d_Jc_edge_fine;  // series current on fine faces (A/m^2)
      extern std::vector<double> sc1d_joule_fine;   // Joule power density on fine cells (W/m^3)
      extern std::vector<int> sc1d_mat_fine;        // majority material of an occupied fine cell, else -1
      // Extra materials in a mixed fine cell. Empty when the cell is pure or empty,
      // so the tunnel weight falls back to sc1d_mat_fine.
      struct fine_mix_t {
         int mat;
         double count;
      };
      extern std::vector<std::vector<fine_mix_t>> sc1d_mix_fine;
      extern std::vector<double> sc1d_sigma_fine;    // electrical conductivity after the interface bake (S/m)
      extern std::vector<double> sc1d_sigma_bulk_fine; // same prolongation before the analytic sheet
      extern std::vector<double> sc1d_seebeck_fine; // Seebeck S (V/K), prolonged from the coarse cell
      extern std::vector<double> sc1d_beta_seebeck_fine; // Seebeck spin factor, not rewritten by the interface sheet
      extern std::vector<double> sc1d_phot_fine; // superdiffusion polarization, not rewritten by the interface sheet
      extern std::vector<double> sc1d_J_seebeck_edge; // Seebeck charge current on fine faces (A/m^2)
      extern std::vector<double> sc1d_J_restore_edge; // J_field = sigma*E on fine faces (A/m^2)
      extern std::vector<double> sc1d_J_hot_edge; // hot-electron charge current, minus the electron flux (A/m^2)
      extern std::vector<double> sc1d_J_cd_edge; // Thomas-Fermi charge diffusion -De d(ne)/dz (A/m^2)
      extern std::vector<double> sc1d_J_sd_edge; // spin-diffusion charge beta_d*De*d(S·m)/dz (A/m^2)
      extern std::vector<double> sc1d_Js_fine;   // total fine-cell spin current, verbose output
      extern std::vector<double> sc1d_mhat_fine; // fine-cell axis, verbose output. Non-magnetic cells are the zero vector.
      extern std::vector<double> sc1d_Te_fine;   // electron temperature on the fine cells (K), one vertical profile
      extern std::vector<double> sc1d_Tp_fine;   // phonon temperature on the fine cells (K)
      extern std::vector<double> sc1d_dTe_dz_edge; // dTe/dz on fine faces (K/m), from the ltmp cell centres
      extern std::vector<double> sc1d_dTp_dz_edge; // dTp/dz on fine faces (K/m)
      extern std::vector<int> sc1d_atom_ls;       // global column index into sc1d_Sfine, -1 if unused
      extern std::vector<int> sc1d_atom_fine;     // occupied fine-cell index of the atom, -1 if unused
      extern std::vector<int> sc1d_atom_fine_lo;  // inclusive fine-cell range
      extern std::vector<int> sc1d_atom_fine_hi;

      // One bit per fine cell of a column. All three are properties of the column,
      // after every rank has added the atoms it holds in that column.
      // Occupied: at least one atom centre lies in the cell.
      // Connected: the cell lies in the field window of such an atom.
      // Spacing: the cell is outside every field window in that column.
      const unsigned char sc1d_mask_occupied = 1;
      const unsigned char sc1d_mask_connected = 2;
      const unsigned char sc1d_mask_spacing = 4;
      extern std::vector<unsigned char> sc1d_cell_mask;
      extern std::vector<int> sc1d_anchor_left;
      extern std::vector<int> sc1d_anchor_right;
      extern std::vector<double> sc1d_anchor_t;
      // Face between fine cells i and i+1. Index is stack * nf + i.
      extern std::vector<double> sc1d_Rint_face;
      extern std::vector<double> sc1d_natom_fine;
      extern std::vector<double> sc1d_msat_fine;
      extern std::vector<double> sc1d_m_prev_fine;
      extern std::vector<double> sc1d_m0_fine;
      extern std::vector<double> sc1d_ml_fine;
      extern std::vector<double> sc1d_mhat_ref_fine;
      extern std::vector<double> sc1d_dmdt_fine;
      extern std::vector<double> sc1d_sheet_acc_up_fine;
      extern std::vector<double> sc1d_sheet_acc_dn_fine;
      extern std::vector<double> sc1d_sheet_acc_re_fine;
      extern std::vector<double> sc1d_sheet_acc_im_fine;
      extern std::vector<double> sc1d_sheet_acc_d2_fine;
      extern std::vector<double> sc1d_Rint_base_face;
      // Atomic planes that own a field window. z in Angstroms. lo < 0 if unused.
      extern std::vector<double> sc1d_plane_z;
      extern std::vector<int> sc1d_plane_lo;
      extern std::vector<int> sc1d_plane_hi;

      // Parallel channel check. A non-physical polarisation aborts.
      double spin_gamma(const double Bc, const double Bd);
      // Occupy the fine grid, prolong from those cells, and bake the interface.
      void build_fine_material_grid();
      // Sum moments onto occupied fine cells and prolong m and dm/dt.
      // The atom map and the anchor slots are fixed for the run.
      void cache_prolong_maps();
      void prolong_fine_magnetisation(const double dt);
      void refresh_fine_transport_constants();
      // Per-atom Eq. 6 coefficients. The fine-grid rates and the atom windows are fixed.
      void cache_sc1d_torque_atoms();

// Internal step counter for the 1D solvers (increments every LLG step when called)
extern unsigned long sc1d_step_counter;

// 0: no chrono, no timing Allgather, no "sc1d timing" log lines.
#ifndef ST_SC1D_TIMINGS
#define ST_SC1D_TIMINGS 0
#endif

// Cumulative wall time (seconds) on this rank for the 1D spin-current step.
extern double sc1d_time_total;
extern double sc1d_time_io;
extern double sc1d_time_charge;
extern double sc1d_time_spin_current;
extern double sc1d_time_spin_acc;
extern double sc1d_time_diffusion;
extern double sc1d_time_eom;
extern double sc1d_time_interp;
extern double sc1d_time_thermal;
extern double sc1d_time_bcast;
extern double sc1d_time_coarsen;
extern double sc1d_time_other;


      // State for demag-driven accumulation and stable-axis handling
      extern std::vector<double> sc1d_m_prev_mag; // |m|(t-dt) per coarse microcell
      extern std::vector<double> sc1d_m0_mag;     // initial |m| per coarse microcell
      extern std::vector<double> sc1d_mhat_ref;   // reference axis per coarse microcell (3*ncells)

      // Fine-grid state per *local* stack (rank-owned stacks only)
      extern int sc1d_nsub;  // fine subdivisions per coarse microcell along z
      extern int sc1d_nf;    // fine nodes per stack == nsub*num_microcells_per_stack
      extern std::vector<int> sc1d_local_stacks;       // list of global stack IDs owned by this rank
      extern std::vector<int> sc1d_stack_local_index;  // size num_stacks_y, maps global stack -> local index or -1
      extern std::vector<double> sc1d_Sfine;           // 3*sc1d_nf per column, every column, on every rank
      extern std::vector<double> sc1d_k_prev;          // Previous RHS for AB2: 3*sc1d_nf per local stack
      // Unused by the one-step TR-BDF2 half-step. Kept so older callers still link.
      extern std::vector<double> sc1d_S_bdf2_h1;
      extern std::vector<double> sc1d_S_bdf2_h2;


      // Persistent fine-grid constant material properties per local stack
      extern std::vector<double> sc1d_Bc_fine;        // beta_cond interpolated to fine grid
      extern std::vector<double> sc1d_Bd_fine;        // beta_diff interpolated to fine grid
      // After the sheet, both equal D_bulk * (σ↑+σ↓) / σ_bulk. beta_c is the imbalance.
      extern std::vector<double> sc1d_D_fine;         // spin and charge diffusivity
      extern std::vector<double> sc1d_Dmat_fine;      // material diffusivity before the sheet; hot-electron speed only
      extern std::vector<double> sc1d_lsf_fine;       // lambda_sdl interpolated to fine grid
      extern std::vector<double> sc1d_lphi_fine;      // lambda_phi interpolated to fine grid
      extern std::vector<double> sc1d_Jsd_fine;       // sd_exchange interpolated to fine grid
      extern std::vector<double> sc1d_chi_fine;       // chi_demag interpolated to fine grid
      extern std::vector<double> sc1d_sa_inf_fine;   // sa_infinity interpolated to fine grid
      extern std::vector<double> sc1d_PN_fine;       // DOS polarization P_N interpolated to fine grid
      extern std::vector<double> sc1d_alpha_fine;    // DOS alpha interpolated to fine grid

      extern std::vector<double> coeff_ast;
      extern std::vector<double> coeff_nast;     
      extern std::vector<double> cell_natom;
      extern std::vector<double> cell_msat;   // sum of atomic moments in the cell, mu_B

      // three-vector arrays
      extern std::vector<double> pos; /// stack position
      extern std::vector<double> m; // magnetisation
      extern std::vector<double> j_final_up_y;
      extern std::vector<double> j_final_up_x;
      extern std::vector<double> j_final_down_y;
      extern std::vector<double> j_int_up_y; // spin current
      extern std::vector<double> j_int_down_y;
      extern std::vector<double> j_init_up_y; // spin current
      extern std::vector<double> j_init_down_y;
      extern std::vector<double> sa_final; // spin accumulation
      extern std::vector<double> ns_final; // non-equilibrium charge density
      extern std::vector<double> ne_final; // excess charge density
      extern std::vector<double> jc_final; // charge current density
      extern std::vector<double> js_final; // spin current density
      // extern std::vector<double> sa_sot_final;
      extern std::vector<double> sa_int;
      // extern std::vector<double> sa_sot_init;
      extern std::vector<double> spin_torque; // spin torque
      extern std::vector<double> ast; // adiabatic spin torque
      extern std::vector<double> nast; // non-adiabatic spin torque
      extern std::vector<double> total_ST; // non-adiabatic spin torque
      extern std::vector<double> magx_mat; // magnetisation of material
      extern std::vector<double> magy_mat;
      extern std::vector<double> magz_mat;
      
      extern std::vector<int> mpi_stack_list_y;
      extern std::vector<int> mpi_stack_list_x;

      //mpi sum variables
      extern std::vector<double> sa_sum;
      extern std::vector<double> ns_sum;
      extern std::vector<double> ne_sum;
      extern std::vector<double> jc_sum;
      extern std::vector<double> js_sum;
      extern std::vector<double> j_final_up_x_sum;
      extern std::vector<double> j_final_up_y_sum;
      extern std::vector<double> j_final_down_y_sum;
      extern std::vector<double> coeff_ast_sum;
      extern std::vector<double> coeff_nast_sum;
      extern std::vector<double> ast_sum;
      extern std::vector<double> nast_sum;
      extern std::vector<double> total_ST_sum;
      extern std::vector<int> cell_natom_sum;


      // material parameters for spin torque calculation
      struct mp_t{
         //STT
         double beta_cond;    /// spin polarisation (conductivity)
         double beta_diff;    /// spin polarisation (diffusion)
         double sa_infinity;  /// intrinsic spin accumulation
         double lambda_sdl;   /// spin diffusion length
         double diffusion;    /// diffusion constant
         double sd_exchange;  /// sd_exchange constant

         // Extensions for 1D transient spin accumulation
         double lambda_phi;   /// transverse spin dephasing length (m)
         double chi_demag;    /// longitudinal demag source on S, (C/m^3) per muB
         double seebeck_coefficient; /// Seebeck coefficient S (V/K) for charge Seebeck effect
         double seebeck_polarization; /// P_S = (S_up - S_down) / S. Default 0.
         double beta_seebeck; /// optional override of the init Seebeck spin factor
         bool beta_seebeck_set; /// true when beta-seebeck was given in the material file
         double superdiffusion_polarization; /// hot-electron polarization. Default 0.
         double conductivity; /// electrical conductivity σ (S/m); <0 unset (Einstein from D, Te); 0 insulator
         double spin_dos_polarization; /// P_N = (N_up - N_down)/(N_up + N_down). Default 0.
         double spin_dos_alpha;        /// alpha = (N_up + N_down)/(2 N_up N_down). Default 1.
         bool spin_currents_thermal;   /// opt into spin-currents-1d-thermal-effects laws

         //SOT
         double sot_beta_cond;    /// spin polarisation (conductivity)
         double sot_beta_diff;    /// spin polarisation (diffusion)
         double sot_sa_infinity;  /// intrinsic spin accumulation
         double sot_lambda_sdl;   /// spin diffusion length
         double sot_diffusion;    /// diffusion constant
         double sot_sd_exchange;  /// sd_exchange constant

         // Thermal gradient properties
         double electron_heat_capacity;        /// electron heat capacity (J/(K^2*m^3))
         double phonon_heat_capacity;          /// phonon heat capacity (J/(K*m^3))
         double electron_thermal_conductivity; /// electron thermal conductivity (W/(m*K))
         double phonon_thermal_conductivity;   /// phonon thermal conductivity (W/(m*K))
         double electron_phonon_coupling;       /// electron-phonon coupling (W/(m^3*K))
         double einstein_temperature;           /// Einstein temperature (K)
      };

      // three vector type definition
      class three_vector_t{
      public:
         double x;
         double y;
         double z;

         // constructor
         three_vector_t(double ix, double iy, double iz){
            x = ix;
            y = iy;
            z = iz;
         }

      };

      // matrix type definition
      struct matrix_t{
         double xx;
         double xy;
         double xz;
         double yx;
         double yy;
         double yz;
         double zx;
         double zy;
         double zz;
      };

      // array of material properties
      extern std::vector<st::internal::mp_t> mp;

      // default material properties
      extern st::internal::mp_t default_properties;

      // -----------------------------------------------------------------------------
      // Optional interfacial spin conductance (Robin-type coupling) support.
      //
      // Users may specify a *pair-wise* interface conductance G_int between material
      // types. Units: m/s. Internally we store the corresponding resistance
      // R_int = 1/G_int with units s/m.
      //
      // - r_int_pair is a dense matrix stored row-major with shape (nmat, nmat)
      //   where nmat == st::internal::mp.size().
      // - r_int_edge is a per-microcell z-edge resistance (coarse grid), stored for
      //   the edge between cell and cell+1 when advancing along z within a stack.
      //   r_int_edge[cell] is valid only when cell is not the last z-layer of a stack.
      //
      // Unspecified pairs default to R_int=0 (i.e. G_int = infinity -> no interface
      // resistance / continuity).
      extern bool interface_coupling_enabled;
      extern std::vector<double> r_int_pair; // size nmat*nmat, units s/m
      extern std::vector<double> r_int_edge; // size ncells, units s/m

      // Ensure r_int_pair has size nmat*nmat (preserving existing values)
      void ensure_interface_matrix_size(const std::size_t nmat);

      // Optional mixing-conductance sheet. Unset pairs stay transparent and
      // leave r_int_edge alone. A set pair replaces that pair's lumped
      // resistance with G↑, G↓, Re/Im G↑↓, and δ on the same coarse edge.
      extern bool sheet_enabled;
      extern std::vector<double> sheet_g_up_pair;    // S/m^2
      extern std::vector<double> sheet_g_down_pair;  // S/m^2
      extern std::vector<double> sheet_re_pair;      // S/m^2
      extern std::vector<double> sheet_im_pair;      // S/m^2
      extern std::vector<double> sheet_delta_pair;   // dimensionless
      extern std::vector<char> sheet_pair_set;       // 1 if that directed pair was parsed
      extern std::vector<double> sheet_g_up_edge;
      extern std::vector<double> sheet_g_down_edge;
      extern std::vector<double> sheet_re_edge;
      extern std::vector<double> sheet_im_edge;
      extern std::vector<double> sheet_delta_edge;
      extern std::vector<char> sheet_edge_on;
      void ensure_sheet_matrix_size(const std::size_t nmat);

      // Tunnel pair, separate from the metallic sheet. Conductances are the
      // zero-thickness channel values (S/m^2). Decay lengths are in metres.
      // A non-positive decay length leaves that channel independent of thickness.
      extern bool tunnel_enabled;
      extern std::vector<double> tunnel_g_up_pair;
      extern std::vector<double> tunnel_g_down_pair;
      extern std::vector<double> tunnel_l_up_pair;
      extern std::vector<double> tunnel_l_down_pair;
      // Unpolarised parallel channel. Prefactors are the 300 K conductances.
      // One follows Te, one follows Tp. A longer decay than the elastic channels
      // makes the inelastic share grow with barrier thickness.
      extern std::vector<double> tunnel_gi_e_pair;
      extern std::vector<double> tunnel_gi_p_pair;
      extern std::vector<double> tunnel_li_pair;
      extern std::vector<double> tunnel_pe_pair;
      extern std::vector<double> tunnel_pp_pair;
      // Bias at which the zero-bias TMR contribution is halved (V).
      // A non-positive value leaves the tunnel pair voltage independent.
      extern std::vector<double> tunnel_vhalf_pair;
      extern std::vector<char> tunnel_pair_set;
      void ensure_tunnel_matrix_size(const std::size_t nmat);

      //-----------------------------------------------------------------------------
      // Shared functions used for the spin torque calculation
      //-----------------------------------------------------------------------------
      void output_microcell_data();
      void output_microcell_sa_data();
      void output_base_microcell_data();
      void output_sc1d_data();
      void output_sc1d_stage_origin(const std::string& stage);
      void output_sc1d_relax_origin();
      void output_sc1d_fine_constants();
      void report_sc1d_timing();
      void report_sc1d_profile(double run_s);
      void calculate_spin_accumulation();
      void calculate_sot_accumulation();
      void initialise_spincurrents_1d();
      void calculate_spin_accumulation_1d();
      void update_cell_magnetisation(const std::vector<double>& x_spin_array,
                                     const std::vector<double>& y_spin_array,
                                     const std::vector<double>& z_spin_array,
                                     const std::vector<int>& atom_type_array,
                                     const std::vector<double>& mu_s_array);

      void set_inverse_transformation_matrix(const st::internal::three_vector_t& reference_vector, st::internal::matrix_t& itm);
      st::internal::three_vector_t transform_vector(const st::internal::three_vector_t& rv, const st::internal::matrix_t& tm);
      st::internal::three_vector_t gaussian_elimination(st::internal::matrix_t& M, st::internal::three_vector_t& V);


      //-----------------------------------------------------------------------------
      // Thermal gradient functions
      //-----------------------------------------------------------------------------
      void initialise_thermal_gradients();
      void update_thermal_gradients_stack(int stack_idx, double time_s, double dt_si);
      void get_thermal_fields(std::vector<double>& thermal_x,
                              std::vector<double>& thermal_y,
                              std::vector<double>& thermal_z,
                              int start_idx, int end_idx);
      
      // Thermal gradient output functions
      void output_thermal_microcell_data();
      void open_thermal_temperature_profile_file();
      void write_thermal_temperature_data();

   } // end of iternal namespace
    extern double spin_acc_time;
} // end of st namespace

#endif //SPINTORQUE_INTERNAL_H_
