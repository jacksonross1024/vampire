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
      extern double sc1d_fine_dz;      // fine grid spacing in Angstroms (along st::internal::stz)
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

// Temperature coupling parameters
extern bool sc1d_use_ltmp_temperatures;        // legacy input; temperature source is ltmp or the global TTM
extern double sc1d_reference_temperature;      // reference temperature for scaling (K)
extern bool sc1d_thermal_effects;              // transport constants may depend on Te, Tp (laws not filled in)

      // Seebeck effect parameters
      extern bool sc1d_seebeck_enable;               // enable charge Seebeck effect

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
      extern std::vector<double> sc1d_ne_seebeck_fine; // Seebeck excess ne (C/m^3), not the Poisson charge
      extern std::vector<double> sc1d_V_fine;        // electrostatic potential (V), size: local_stacks*nf
      extern std::vector<double> sc1d_Jc_edge_fine;  // Jc faces (A/m^2), size: local_stacks*(nf+1)
      extern std::vector<double> sc1d_sigma_fine;    // electrical conductivity (S/m), prolonged from the coarse cell
      extern std::vector<double> sc1d_seebeck_fine;  // Seebeck S (V/K), prolonged from the coarse cell
      extern std::vector<int>    sc1d_mat_fine;      // unused; transport no longer picks a material per fine cell
      extern std::vector<double> sc1d_sa_demag_fine; // legacy storage; target channel removed

      // Atom ↔ fine-grid map for conservative S → atom torque
      extern std::vector<int> sc1d_atom_ls;       // global column index into sc1d_Sfine, -1 if unused
      extern std::vector<int> sc1d_atom_fine_lo;  // inclusive fine-cell range
      extern std::vector<int> sc1d_atom_fine_hi;

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
extern double sc1d_time_interp;
extern double sc1d_time_bcast;
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


      // Persistent fine-grid constant material properties per local stack
      extern std::vector<double> sc1d_Bc_fine;        // beta_cond interpolated to fine grid
      extern std::vector<double> sc1d_Bd_fine;        // beta_diff interpolated to fine grid
      extern std::vector<double> sc1d_D_fine;         // diffusion interpolated to fine grid
      extern std::vector<double> sc1d_lsf_fine;       // lambda_sdl interpolated to fine grid
      extern std::vector<double> sc1d_lphi_fine;      // lambda_phi interpolated to fine grid
      extern std::vector<double> sc1d_Jsd_fine;       // sd_exchange interpolated to fine grid
      extern std::vector<double> sc1d_chi_fine;       // chi_demag interpolated to fine grid
      extern std::vector<double> sc1d_sa_inf_fine;   // sa_infinity interpolated to fine grid
      extern std::vector<double> sc1d_alpha_edge;     // diffusion operator edge coefficients

      extern std::vector<double> coeff_ast;
      extern std::vector<double> coeff_nast;     
      extern std::vector<double> cell_natom;

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
         double conductivity; /// electrical conductivity σ (S/m); <0 unset (Einstein from D, Te); 0 insulator

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

      //-----------------------------------------------------------------------------
      // Shared functions used for the spin torque calculation
      //-----------------------------------------------------------------------------
      void output_microcell_data();
      void output_microcell_sa_data();
      void output_base_microcell_data();
      void output_sc1d_data();
      void report_sc1d_timing();
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
