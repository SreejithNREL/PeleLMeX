#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include <PeleLMDeriveFunc.H>
#include <EOS.H>
#include <Transport.H>

using namespace amrex;

static Box the_same_box (const Box& b)    { return b;                 }

void PeleLM::Setup() {
   BL_PROFILE_VAR("PeleLM::Setup()", Setup);

   // Read PeleLM parameters
   readParameters();

   // Setup the state variables
   variablesSetup();

   // Derived variables
   derivedSetup();

   // Tagging setup
   taggingSetup();

   // Initialize Level Hierarchy data
   resizeArray();

   // Initiliaze BCs
   setBoundaryConditions();

   // Initialize EOS and others
   if (!m_incompressible) {
      amrex::Print() << " Initialization of EOS ... \n";
      EOS::init();
      amrex::Print() << " Initialization of Transport ... \n";
      transport_init();
   }

   // Problem parameters
   prob_parm.reset(new ProbParm{});

   // Problem parameters
   readProbParm();
}

void PeleLM::readParameters() {
   BL_PROFILE_VAR("PeleLM::readParameters()", readParameters);

   readIOParameters();

   ParmParse pp("peleLM");

   // -----------------------------------------
   // Misc
   // -----------------------------------------
   pp.query("v", m_verbose);

   // -----------------------------------------
   // Boundary conditions
   // -----------------------------------------

   Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
   Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
   pp.getarr("lo_bc",lo_bc_char,0,AMREX_SPACEDIM);
   pp.getarr("hi_bc",hi_bc_char,0,AMREX_SPACEDIM);

   Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
   for (int dir = 0; dir<AMREX_SPACEDIM; dir++){
      if (lo_bc_char[dir] == "Interior") {
         lo_bc[dir] = 0;
      } else if (lo_bc_char[dir] == "Inflow") {
         lo_bc[dir] = 1;
      } else if (lo_bc_char[dir] == "Outflow") {
         lo_bc[dir] = 2;
      } else if (lo_bc_char[dir] == "Symmetry") {
         lo_bc[dir] = 3;
      } else if (lo_bc_char[dir] == "SlipWallAdiab") {
         lo_bc[dir] = 4;
      } else if (lo_bc_char[dir] == "NoSlipWallAdiab") {
         lo_bc[dir] = 5;
      } else if (lo_bc_char[dir] == "SlipWallIsotherm") {
         lo_bc[dir] = 6;
      } else if (lo_bc_char[dir] == "NoSlipWallIsotherm") {
         lo_bc[dir] = 7;
      } else {
         amrex::Abort("Wrong boundary condition word in lo_bc, please use: Interior, Inflow, Outflow, "
                      "Symmetry, SlipWallAdiab, NoSlipWallAdiab, SlipWallIsotherm, NoSlipWallIsotherm");
      }

      if (hi_bc_char[dir] == "Interior") {
         hi_bc[dir] = 0;
      } else if (hi_bc_char[dir] == "Inflow") {
         hi_bc[dir] = 1;
      } else if (hi_bc_char[dir] == "Outflow") {
         hi_bc[dir] = 2;
      } else if (hi_bc_char[dir] == "Symmetry") {
         hi_bc[dir] = 3;
      } else if (hi_bc_char[dir] == "SlipWallAdiab") {
         hi_bc[dir] = 4;
      } else if (hi_bc_char[dir] == "NoSlipWallAdiab") {
         hi_bc[dir] = 5;
      } else if (hi_bc_char[dir] == "SlipWallIsotherm") {
         hi_bc[dir] = 6;
      } else if (hi_bc_char[dir] == "NoSlipWallIsotherm") {
         hi_bc[dir] = 7;
      } else {
         amrex::Abort("Wrong boundary condition word in hi_bc, please use: Interior, Inflow, Outflow, "
                      "Symmetry, SlipWallAdiab, NoSlipWallAdiab, SlipWallIsotherm, NoSlipWallIsotherm");
      }
   }

   // Store BCs in m_phys_bc.
   for (int i = 0; i < AMREX_SPACEDIM; i++) {
      m_phys_bc.setLo(i,lo_bc[i]);
      m_phys_bc.setHi(i,hi_bc[i]);
   }

   // -----------------------------------------
   // Algorithm
   // -----------------------------------------
   pp.query("use_divu", m_has_divu);
   pp.query("incompressible", m_incompressible);
   if (m_incompressible) m_has_divu = 0;
   pp.query("rho", m_rho);
   pp.query("mu", m_mu);
}

void PeleLM::readIOParameters() {
   BL_PROFILE_VAR("PeleLM::readIOParameters()", readIOParameters);

   ParmParse pp("amr");

   pp.query("plot_file", m_plot_file);
   pp.query("plot_int" , m_plot_int);
   m_derivePlotVarCount = (pp.countval("derive_plot_vars"));
   if (m_derivePlotVarCount != 0) {
      m_derivePlotVars.resize(m_derivePlotVarCount);
      for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++) {
         pp.get("derive_plot_vars", m_derivePlotVars[ivar],ivar);
      }
   }

   pp.query("max_step", m_max_step);
   pp.query("stop_time", m_stop_time);
   pp.query("fixed_dt", m_fixed_dt);
   pp.query("cfl", m_cfl);
   pp.query("dt_shrink", m_dtshrink);
   pp.query("dt_change_max", m_dtChangeMax);

}

void PeleLM::variablesSetup() {
   BL_PROFILE_VAR("PeleLM::variablesSetup()", variablesSetup);

   std::string PrettyLine = std::string(78, '=') + "\n";

   //----------------------------------------------------------------
   // Variables ordering is defined through compiler macro in PeleLM_Index.H 
   // Simply print on screen the state layout and append to the stateComponents list
   Print() << "\n";
   Print() << PrettyLine;
   Print() << " State components \n";
   Print() << PrettyLine;

   Print() << " Velocity X: " << VELX;
   stateComponents.emplace_back(VELX,"x_velocity");
#if AMREX_SPACEDIM > 1
   Print() << ", Velocity Y: " << VELY;
   stateComponents.emplace_back(VELY,"y_velocity");
#if AMREX_SPACEDIM > 2
   Print() << ", Velocity Z: " << VELZ;
   stateComponents.emplace_back(VELZ,"z_velocity");
#endif
#endif
   Print() << " \n";

   if (! m_incompressible) {
      Print() << " Density: " << DENSITY << "\n";
      stateComponents.emplace_back(DENSITY,"density");
      Print() << " First species: " << FIRSTSPEC << "\n";
      Vector<std::string> names;
      EOS::speciesNames(names);
      for (int n = 0; n < NUM_SPECIES; n++ ) {
         stateComponents.emplace_back(FIRSTSPEC+n,"rho.Y("+names[n]+")");
      }
      Print() << " Enthalpy: " << RHOH << "\n";
      stateComponents.emplace_back(RHOH,"rhoh");
      Print() << " Temperature: " << TEMP << "\n";
      stateComponents.emplace_back(TEMP,"temp");
      Print() << " thermo. pressure: " << RHORT << "\n";
      stateComponents.emplace_back(RHORT,"RhoRT");
   }

   if (m_nAux > 0) {
      Print() << " First passive scalar: " << FIRSTAUX << "\n";
      for (int n = 0; n < m_nAux; n++ ) {
         stateComponents.emplace_back(FIRSTAUX+n,"Aux_"+std::to_string(n));
      }
   }

   Print() << " => Total number of state variables: " << NVAR << "\n";
   Print() << PrettyLine;
   Print() << "\n";

   //----------------------------------------------------------------
   // Set advection/diffusion types
   m_AdvTypeState.resize(NVAR);
   m_DiffTypeState.resize(NVAR);

   // Velocity - follow incflo
   for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
#ifdef AMREX_USE_EB
      m_AdvTypeState[VELX+idim] = 1;   // Conservative
#else
      m_AdvTypeState[VELX+idim] = 0;   // NonConservative
#endif
      m_DiffTypeState[VELX+idim] = 1;  // Diffusive
   }

   if (! m_incompressible) {
      m_AdvTypeState[DENSITY] = 1;
      m_DiffTypeState[DENSITY] = 0;
      for (int n = 0; n < NUM_SPECIES; ++n) {
         m_AdvTypeState[FIRSTSPEC+n] = 1;
         m_DiffTypeState[FIRSTSPEC+n] = 1;
      }
      m_AdvTypeState[RHOH] = 1;
      m_DiffTypeState[RHOH] = 1;
      m_AdvTypeState[TEMP] = 0;
      m_DiffTypeState[TEMP] = 0;
      m_AdvTypeState[RHORT] = 0;
      m_DiffTypeState[RHORT] = 0;
   }
}

void PeleLM::derivedSetup()
{
   BL_PROFILE_VAR("PeleLM::derivedSetup()", derivedSetup);

   if (!m_incompressible) {

      // Get species names
      Vector<std::string> spec_names;
      EOS::speciesNames(spec_names);

      // Set species mass fractions
      Vector<std::string> var_names_massfrac(NUM_SPECIES);
      for (int n = 0 ; n < NUM_SPECIES; n++) {
         var_names_massfrac[n] = "Y("+spec_names[n]+")";
      }
      derive_lst.add("mass_fractions",IndexType::TheCellType(),NUM_SPECIES,
                     var_names_massfrac,pelelm_dermassfrac,the_same_box);

   }

   // Cell average pressure
   derive_lst.add("avg_pressure",IndexType::TheCellType(),1,pelelm_deravgpress,the_same_box);

}

void PeleLM::taggingSetup()
{
   BL_PROFILE_VAR("PeleLM::taggingSetup()", taggingSetup);

   std::string amr_prefix = "amr";
   ParmParse ppamr(amr_prefix);

   Vector<std::string> refinement_indicators;
   ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));
   for (int n = 0; n<refinement_indicators.size(); ++n) {
      std::string ref_prefix = amr_prefix + "." + refinement_indicators[n];
      ParmParse ppr(ref_prefix);

      // Tag a given box
      RealBox realbox;
      if (ppr.countval("in_box_lo")) {
         Vector<Real> box_lo(AMREX_SPACEDIM);
         Vector<Real> box_hi(AMREX_SPACEDIM);
         ppr.getarr("in_box_lo",box_lo,0,box_lo.size());
         ppr.getarr("in_box_hi",box_hi,0,box_hi.size());
         realbox = RealBox(&(box_lo[0]),&(box_hi[0]));
      }

      AMRErrorTagInfo info;

      if (realbox.ok()) {
         info.SetRealBox(realbox);
      }

      if (ppr.countval("start_time") > 0) {
         Real min_time; ppr.get("start_time",min_time);
         info.SetMinTime(min_time);
      }

      if (ppr.countval("end_time") > 0) {
         Real max_time; ppr.get("end_time",max_time);
         info.SetMaxTime(max_time);
      }

      if (ppr.countval("max_level") > 0) {
         int max_level; ppr.get("max_level",max_level);
         info.SetMaxLevel(max_level);
      }

      bool itexists = false;
      if (ppr.countval("value_greater")) {
         Real value; ppr.get("value_greater",value);
         std::string field; ppr.get("field_name",field);
         errTags.push_back(AMRErrorTag(value,AMRErrorTag::GREATER,field,info));
         itexists = derive_lst.canDerive(field) || isStateVariable(field);
      } else if (ppr.countval("value_less")) {
         Real value; ppr.get("value_less",value);
         std::string field; ppr.get("field_name",field);
         errTags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,field,info));
         itexists = derive_lst.canDerive(field) || isStateVariable(field);
      } else if (ppr.countval("vorticity_greater")) {
         Real value; ppr.get("vorticity_greater",value);
         const std::string field="mag_vort";
         errTags.push_back(AMRErrorTag(value,AMRErrorTag::VORT,field,info));
         itexists = derive_lst.canDerive(field) || isStateVariable(field);
      } else if (ppr.countval("adjacent_difference_greater")) {
         Real value; ppr.get("adjacent_difference_greater",value);
         std::string field; ppr.get("field_name",field);
         errTags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,field,info));
         itexists = derive_lst.canDerive(field) || isStateVariable(field);
      } else if (realbox.ok()) {
        errTags.push_back(AMRErrorTag(info));
      } else {
        Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[n]).c_str());
      }

      if ( !itexists ) {
         amrex::Error("PeleLM::taggingSetup(): unknown variable field for criteria "+refinement_indicators[n]);
      }
   }
}

void PeleLM::resizeArray() {

   if (m_verbose) {
      Print() << " Initializing data for " << max_level+1 << " levels \n";
   }

   // State data
   m_leveldata_old.resize(max_level+1);
   m_leveldata_new.resize(max_level+1);

   // Factory
   m_factory.resize(max_level+1);

   // Time
   m_t_old.resize(max_level+1);
   m_t_new.resize(max_level+1);

}
