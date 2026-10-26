#include <PeleLMeX.H>
#include <algorithm>
#include <map>

void
PeleLM::createDiagnostics()
{
  std::string pele_prefix = "peleLM";
  amrex::ParmParse pp(pele_prefix);
  int n_diags = 0;
  n_diags = pp.countval("diagnostics");
  amrex::Vector<std::string> diags;
  if (n_diags > 0) {
    m_diagnostics.resize(n_diags);
    diags.resize(n_diags);
  }
  for (int n = 0; n < n_diags; ++n) {
    pp.get("diagnostics", diags[n], n);
    std::string diag_prefix = pele_prefix + "." + diags[n];
    amrex::ParmParse ppd(diag_prefix);
    std::string diag_type;
    ppd.get("type", diag_type);
    m_diagnostics[n] = DiagBase::create(diag_type);
    m_diagnostics[n]->init(diag_prefix, m_base_output_pref + diags[n]);
    m_diagnostics[n]->addVars(m_diagVars);
  }

  // Remove duplicates from m_diagVars and check that all the variables exists
  std::sort(m_diagVars.begin(), m_diagVars.end());
  auto last = std::unique(m_diagVars.begin(), m_diagVars.end());
  m_diagVars.erase(last, m_diagVars.end());
  for (auto& v : m_diagVars) {
    bool itexists =
      derive_lst.canDerive(v) || isStateVariable(v) || isReactVariable(v);
    if (!itexists) {
      amrex::Abort("Field " + v + " is not available");
    } else {
      if (derive_lst.canDerive(v)) {
        const PeleLMDeriveRec* rec = derive_lst.get(v);
        if (rec->variableComp(v) < 0) {
          std::string errmsg = "Diagnostics can't handle derived with more "
                               "than 1 component at the moment.\n";
          errmsg += "Add the desired components individually.\n";
          amrex::Abort(errmsg);
        }
      }
    }
  }
}

void
PeleLM::updateDiagnostics()
{
  // Might need to update some internal data as the grid changes
  for (const auto& m_diagnostic : m_diagnostics) {
    if (m_diagnostic->needUpdate()) {
      m_diagnostic->prepare(
        finestLevel() + 1, Geom(0, finestLevel()), boxArray(0, finestLevel()),
        dmap, m_diagVars);
    }
  }
}

void
PeleLM::doDiagnostics()
{
  BL_PROFILE("PeleLMeX::doDiagnostics()");

  // Only assemble the variables needed by the diagnostics due this step
  amrex::Vector<int> dueDiags;
  amrex::Vector<std::string> neededVars;
  for (int n = 0; n < static_cast<int>(m_diagnostics.size()); ++n) {
    if (m_diagnostics[n]->doDiag(m_cur_time, m_nstep)) {
      dueDiags.push_back(n);
      m_diagnostics[n]->addVars(neededVars);
    }
  }
  if (dueDiags.empty()) {
    return;
  }
  std::sort(neededVars.begin(), neededVars.end());
  neededVars.erase(
    std::unique(neededVars.begin(), neededVars.end()), neededVars.end());

  // Assemble a vector of MF containing the requested data, in the order of
  // m_diagVars (the diagnostics index into it); unneeded components stay 0.
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> diagMFVec(finestLevel() + 1);
  for (int lev{0}; lev <= finestLevel(); ++lev) {
    diagMFVec[lev] = std::make_unique<amrex::MultiFab>(
      grids[lev], dmap[lev], m_diagVars.size(), 1);
    diagMFVec[lev]->setVal(0.0);
    // Each derive record is evaluated at most once per level, even when
    // several of its components are requested
    std::map<const PeleLMDeriveRec*, std::unique_ptr<amrex::MultiFab>> derCache;
    for (int v{0}; v < m_diagVars.size(); ++v) {
      if (!std::binary_search(
            neededVars.begin(), neededVars.end(), m_diagVars[v])) {
        continue;
      }
      const PeleLMDeriveRec* rec = derive_lst.get(m_diagVars[v]);
      if (rec != nullptr) {
        auto it = derCache.find(rec);
        if (it == derCache.end()) {
          it = derCache.emplace(rec, derive(rec->name(), m_cur_time, lev, 1))
                 .first;
        }
        const int mf_idx = rec->variableComp(m_diagVars[v]);
        amrex::MultiFab::Copy(*diagMFVec[lev], *(it->second), mf_idx, v, 1, 1);
      } else {
        std::unique_ptr<amrex::MultiFab> mf =
          derive(m_diagVars[v], m_cur_time, lev, 1);
        amrex::MultiFab::Copy(*diagMFVec[lev], *mf, 0, v, 1, 1);
      }
    }
  }

  for (const int n : dueDiags) {
    m_diagnostics[n]->processDiag(
      m_nstep, m_cur_time, GetVecOfConstPtrs(diagMFVec), m_diagVars);
  }
}
