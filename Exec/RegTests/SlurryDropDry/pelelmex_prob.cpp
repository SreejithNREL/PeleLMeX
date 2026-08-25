#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

namespace {
// Antoine equation for water, NIST coefficients valid over 344 - 373 K
// (extrapolated slightly to 374.15 K).  log10(p[bar]) = A - B/(T + C)
amrex::Real
psat_water(const amrex::Real T)
{
  constexpr amrex::Real A = 5.0768;
  constexpr amrex::Real B = 1659.793;
  constexpr amrex::Real C = -45.854;
  return 1.e5 * std::pow(10., A - B / (T + C)); // Pa
}
} // namespace

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  PeleLM::prob_parm->eosparm = PeleLM::eos_parms.device_parm();

  // Gas phase properties
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("T0_gas", PeleLM::prob_parm->T0_gas);
  pp.query("vel_gas", PeleLM::prob_parm->vel_gas);
  pp.query("N2_gas", PeleLM::prob_parm->Y_N2);
  pp.query("O2_gas", PeleLM::prob_parm->Y_O2);
  pp.query("RH", PeleLM::prob_parm->RH);

  // Particle properties (only used for reporting here; the spray container
  // reads these itself in SprayParticlesInitInsert.cpp)
  amrex::Real drop_dia = 0.;
  amrex::Real T_d = 0.;
  pp.query("part_dia", drop_dia);
  pp.query("part_temp", T_d);

  // -----------------------------------------------------------------------
  // Convert relative humidity into a water vapor mass fraction and rescale
  // the dry air mass fractions so the mixture sums to one.
  // -----------------------------------------------------------------------
  const amrex::Real T_g = PeleLM::prob_parm->T0_gas;
  const amrex::Real p_g = PeleLM::prob_parm->P_mean;
  const amrex::Real RH = PeleLM::prob_parm->RH;
  amrex::Real Y_H2O = 0.;
  if (RH > 0.) {
    const amrex::Real p_v = RH * psat_water(T_g);
    const amrex::Real X_v = amrex::min<amrex::Real>(1., p_v / p_g);
    // Dry air molecular weight
    constexpr amrex::Real mw_air = 28.85e-3;  // kg/mol
    constexpr amrex::Real mw_h2o = 18.015e-3; // kg/mol
    Y_H2O = X_v * mw_h2o / (X_v * mw_h2o + (1. - X_v) * mw_air);
  }
  PeleLM::prob_parm->Y_H2O = Y_H2O;
  // Rescale dry air composition
  {
    const amrex::Real sum_dry =
      PeleLM::prob_parm->Y_N2 + PeleLM::prob_parm->Y_O2;
    PeleLM::prob_parm->Y_N2 *= (1. - Y_H2O) / sum_dry;
    PeleLM::prob_parm->Y_O2 *= (1. - Y_H2O) / sum_dry;
  }

  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream ofs("ic.txt", std::ofstream::out);
    amrex::Print(ofs) << "T0_gas = " << T_g << "\n"
                      << "P_mean = " << p_g << "\n"
                      << "RH = " << RH << "\n"
                      << "psat_H2O = " << psat_water(T_g) << "\n"
                      << "Y_H2O = " << PeleLM::prob_parm->Y_H2O << "\n"
                      << "Y_O2 = " << PeleLM::prob_parm->Y_O2 << "\n"
                      << "Y_N2 = " << PeleLM::prob_parm->Y_N2 << "\n"
                      << "vel_gas = " << PeleLM::prob_parm->vel_gas << "\n"
                      << "drop_dia = " << drop_dia << "\n"
                      << "drop_temp = " << T_d << std::endl;
    ofs.close();
  }

  // Read mesh-mapping scaling factors for spray IC
  {
    amrex::ParmParse ppg("geometry");
    std::string mesh_map = "ConstantMap";
    if (ppg.countval("mesh_mapping") > 0) {
      ppg.get("mesh_mapping", mesh_map);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        mesh_map == "ConstantMap",
        "Sprays can only be mapped with ConstantMap");
      amrex::ParmParse ppcm("ConstantMap");
      amrex::Vector<amrex::Real> fac(AMREX_SPACEDIM, 1.0);
      ppcm.queryarr("scaling_factor", fac, 0, AMREX_SPACEDIM);
      AMREX_D_TERM(prob_parm->fac_x = fac[0];, prob_parm->fac_y = fac[1];
                   , prob_parm->fac_z = fac[2];);
    }
  }
}

void
PeleLM::freeProbParm()
{
}
