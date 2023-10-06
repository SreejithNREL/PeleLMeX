#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   bool activate_fuel_injection = false;
   amrex::Real fuel_injection_max_radius=0.5;
   //amrex::Real Y_O2_FUEL=0.065963;
   //amrex::Real Y_N2_FUEL=0.434037;
   //amrex::Real Y_C12H26_FUEL=0.5;

   //mass flow correction factors
   amrex::Real mfr_factor_pilot=1.0;
   amrex::Real mfr_factor_cp00=1.0;
   amrex::Real mfr_factor_cp01=1.0;
   amrex::Real mfr_factor_cp02=1.0;
   amrex::Real mfr_factor_cp03=1.0;
   amrex::Real mfr_factor_bch=1.0;
   
   // Chamber conditions
   amrex::Real P_mean        = 101325.0;
   amrex::Real T_mean        = 300.0;

   // Premixers
   amrex::Real premix_phi = 1.136;
   amrex::Real premix_fuel_air_temperature_cp00 = 300;		//These are the fuel-air premixed gas temperature. Since the fuel mass fractions are low, I am assuming it to be T3
   amrex::Real premix_fuel_air_temperature_cp01 = 300;          //If in future, the temperature needs to be altered, these params could be used. All temperatures in Kelvin
   amrex::Real premix_fuel_air_temperature_cp02 = 300;
   amrex::Real premix_fuel_air_temperature_cp03 = 300;
   amrex::Real premix_fuel_air_temperature_pilot = 300;

   // Ignition
   int do_ignit = 0;
   amrex::Real ignit_rad = 0.008;
   amrex::Real ignit_lowZ = 0.02;
   amrex::Real ignit_temp = 1400.0;
   
};
#endif
