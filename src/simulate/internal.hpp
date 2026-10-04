#ifndef SIM_INTERNAL_H_
#define SIM_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// simulation methods implementation. These functions should
// not be accessed outside of the simulate module.
//---------------------------------------------------------------------

namespace sim{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal data types used for simulation module
      //-----------------------------------------------------------------------------

      // simple initialised class for set variables
      class set_double_t{

      private:
         double value; // value
         bool setf; // flag specifiying variable has been set

      public:
         // class functions
         // constructor
         set_double_t() : value(0.0), setf(false) { }

         // setting function
         void set(double in_value){
            value = in_value;
            setf = true;
         };

         // get value function
         double get(){ return value; };
         // check if variable is set
         bool is_set(){ return setf; };

      };

      struct mp_t{
         set_double_t stt_asm; // spin tranfer torque asymmetry
         set_double_t stt_rj;  // spin tranfer relaxation torque
         set_double_t stt_pj;  // spin transfer precession torque
         set_double_t sot_asm; // spin orbit torque asymmetry
         set_double_t sot_asm_2nd_order; // spin orbit torque asymmetry
         set_double_t sot_rj;  // spin orbit relaxation torque
         set_double_t sot_pj;  // spin orbit precession torque
         set_double_t sot_asm2; // spin orbit torque asymmetry
         set_double_t sot_rj2;  // spin orbit relaxation torque
         set_double_t sot_pj2;
         set_double_t let_stt_fl; // laser-electrical STT field-like torque (Tesla at Jc_ref)
         set_double_t let_stt_dl; // laser-electrical STT damping-like torque (Tesla at Jc_ref)
         set_double_t let_fl_torkance; // field-like torkance (Tesla / (A/m^2))
         set_double_t let_dl_torkance; // damping-like torkance (Tesla / (A/m^2))
         set_double_t let_diffusion; // De (m^2/s), shared by field-like and damping-like Eq. 6
         set_double_t let_lambda_J; // exchange rotation length λJ (m)
         set_double_t let_lambda_phi; // spin dephasing length λφ (m)
         set_double_t let_abs_factor; // S built over k λ; default 1
         set_double_t let_P; // Eq. 4 current spin polarisation P
         bool let_pvec_set; // material-level let-polarisation-vector
         double let_pvec[3];
         set_double_t vcmak;   // voltage controlled anisotropy coefficient
         set_double_t lt_x;
         set_double_t lt_y;
         set_double_t lt_z;
         mp_t() : let_pvec_set(false) {
            let_pvec[0] = 0.0;
            let_pvec[1] = 0.0;
            let_pvec[2] = 0.0;
         }
      };

      //-----------------------------------------------------------------------------
      // Internal shared variables used for the simulation
      //-----------------------------------------------------------------------------
      extern bool enable_spin_torque_fields; // flag to enable spin torque fields
      extern bool enable_vcma_fields;        // flag to enable voltage-controlled anisotropy fields

      extern std::vector<sim::internal::mp_t> mp; // array of material properties

      extern std::vector<double> stt_asm; // array of spin transfer torque asymmetry
      extern std::vector<double> stt_rj; // array of adiabatic spin torques
      extern std::vector<double> stt_pj; // array of non-adiabatic spin torques
      extern std::vector<double> stt_polarization_unit_vector; // stt spin polarization direction

      extern std::vector<double> sot_asm; // array of spin orbit torque asymmetry
      extern std::vector<double> sot_asm_2nd_order; // array of spin orbit torque asymmetry
      extern std::vector<double> sot_rj; // array of adiabatic spin torques
      extern std::vector<double> sot_pj; // array of non-adiabatic spin torques
      extern std::vector<double> sot_asm2; // array of spin orbit torque asymmetry
      extern std::vector<double> sot_rj2; // array of adiabatic spin torques
      extern std::vector<double> sot_pj2;
      extern std::vector<double> sot_polarization_unit_vector; // sot spin polarization direction
      extern std::vector<double> sot_polarization_unit_vector2;
      extern std::vector<double> let_stt_fl; // laser-electrical STT field-like torque (Tesla at Jc_ref)
      extern std::vector<double> let_stt_dl; // laser-electrical STT damping-like torque (Tesla at Jc_ref)
      extern std::vector<double> let_fl_torkance; // Tesla / (A/m^2); scaled by Jc
      extern std::vector<double> let_dl_torkance; // Tesla / (A/m^2); scaled by Jc
      extern std::vector<double> let_eq6_fl; // Eq. 6 field-like: Tesla per Js (A/s)
      extern std::vector<double> let_eq6_dl; // Eq. 6 damping-like: Tesla per Js (A/s)
      extern std::vector<int> let_fl_mode; // 0 Tesla*Jc/Jc_ref, 1 torkance*Js, 2 Eq.6*Js
      extern std::vector<int> let_dl_mode;
      extern std::vector<double> let_P; // Eq. 4 current spin polarisation P (per material)
      extern std::vector<double> let_px; // per-material LET polarisation
      extern std::vector<double> let_py;
      extern std::vector<double> let_pz;
      extern double electrical_pulse_strength;

      extern std::vector<double> lot_lt_x;
      extern std::vector<double> lot_lt_y;
      extern std::vector<double> lot_lt_z;

      extern std::vector<double> lot_unit_vector;
      extern std::vector<double> vcmak;   // voltage controlled anisotropy coefficient
      
      // shared Functions
      void llg_quantum_step();

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      extern void initialize_modules();
      extern void increment_time();

   //MPI variables
       extern std::vector<std::vector<int> > c_octants; //Core atoms of each octant
       extern std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant
   } // end of internal namespace
  
} // end of sim namespace

#endif //SIM_INTERNAL_H_
