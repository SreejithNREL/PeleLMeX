#include "PeleLMeX_BPatch.H"
#include "PelePhysics.H"


BPatch::BPatch(const std::string& patch_name,const amrex::Geometry& geom):m_patchname(std::move(patch_name)){

	std::string ppbpatch = "bpatch." + m_patchname;
	amrex::ParmParse ps(ppbpatch);

	const auto dx = geom.CellSizeArray();
	auto prob_lo = geom.ProbLoArray();
	auto prob_hi = geom.ProbHiArray();

	ps.get("patchtype", m_patchtype);

	if(m_patchtype!="fullboundary" && m_patchtype!="circle" && m_patchtype!="rectangle" && m_patchtype!="circle-annular" && m_patchtype!="rectangule-annular"){
		amrex::Abort("\nPatch type do not match allowed values (fullboundary,circle,rectangle,circle-annular,rectangular-annular) \n");
		amrex::Print()<<"\nPatch type = "<<m_patchtype;
	}
	ps.get("boundary_dimension", m_bpdata_h.m_boundary_dim);
	if(m_bpdata_h.m_boundary_dim!=0 && m_bpdata_h.m_boundary_dim!=1 && m_bpdata_h.m_boundary_dim!=2)
	{
		amrex::Abort("\nBoundary dimension should be 0,1 or 2");
	}
	ps.get("low_or_high", m_bpdata_h.m_boundary_lo_hi);
	if(m_bpdata_h.m_boundary_lo_hi!=0 && m_bpdata_h.m_boundary_lo_hi!=1)
	{
		amrex::Abort("\nBoundary high low should be 0 or 1");
	}
	if(m_patchtype=="fullboundary")
	{
		m_bpdata_h.m_patchtype_num=0;
		//treat it like its a rectangular boundary

	}
	else if(m_patchtype=="circle"){
		m_bpdata_h.m_patchtype_num=0;
		ps.get("circle_radius", m_bpdata_h.m_patch_circle_radius);
		if(m_bpdata_h.m_patch_circle_radius<=0.0)
		{
			amrex::Abort("\nPatch type of circle should have radius greater than 0");

		}
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
		    ps.get("circle_center", m_bpdata_h.m_patch_circle_center[n], n);
		}
	}
	else if(m_patchtype=="rectangle"){
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
			ps.get("patch_rectangle_lo", m_bpdata_h.m_patch_rectangle_lo[n], n);
		}
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
			ps.get("patch_rectangle_hi", m_bpdata_h.m_patch_rectangle_hi[n], n);
		}
	}
	else if(m_patchtype=="circle-annular"){
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
			ps.get("patch_circ_ann_center", m_bpdata_h.m_patch_circ_ann_center[n], n);
		}

		ps.get("patch_circ_ann_inner_radius",m_bpdata_h.m_patch_circ_ann_inner_radius);
		ps.get("patch_circ_ann_outer_radius",m_bpdata_h.m_patch_circ_ann_outer_radius);
		if(m_bpdata_h.m_patch_circ_ann_outer_radius<=0.0 || m_bpdata_h.m_patch_circ_ann_inner_radius<=0.0 || m_bpdata_h.m_patch_circ_ann_inner_radius>=m_bpdata_h.m_patch_circ_ann_outer_radius)
		{
			amrex::Abort("\nRadii do not look good for annular circles");
		}
	}
	else if(m_patchtype=="rectangule-annular"){
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
			ps.get("patch_rect_ann_outer_lo", m_bpdata_h.m_patch_rect_ann_outer_lo[n], n);
		}
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
			ps.get("patch_rect_ann_outer_hi", m_bpdata_h.m_patch_rect_ann_outer_hi[n], n);
		}
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
			ps.get("patch_rect_ann_inner_lo", m_bpdata_h.m_patch_rect_ann_inner_lo[n], n);
		}
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
			ps.get("patch_rect_ann_inner_hi", m_bpdata_h.m_patch_rect_ann_inner_hi[n], n);
		}
	}
	else{
		amrex::Abort("\nError! Unknown patch type");
	}

	m_bpdata_h.num_species = ps.countval("species");
	ps.getarr("species", speciesList);

	if (m_bpdata_h.num_species > 0){
		speciesList.resize(m_bpdata_h.num_species);
		m_bpdata_h.speciesIndex = (int*)amrex::The_Pinned_Arena()->alloc(
				m_bpdata_h.num_species * sizeof(int));
		m_bpdata_h.speciesFlux = (amrex::Real*)amrex::The_Pinned_Arena()->alloc(
				m_bpdata_h.num_species * sizeof(amrex::Real));
	}
	else
	{
		amrex::Abort("\nError! No species provided to plot flux");
	}

	amrex::Vector<std::string> names;
	pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(names);
	names.resize(names.size());

	for(int n =0;n<names.size();n++){
		m_bpdata_h.speciesIndex[n]=-1;
	}

	for(int m=0;m<m_bpdata_h.num_species;m++){
		for(int n =0;n<names.size();n++){
			if(speciesList[m]==names[n])
			{
				m_bpdata_h.speciesIndex[m]=n;
			}
		}

	}

	for(int n =0;n<m_bpdata_h.num_species;n++){
		if(m_bpdata_h.speciesIndex[n]==-1){
			std::string msg = "\nError! Unable to find species index "+std::to_string(n);
			amrex::Abort(msg);
		}
	}

	allocate();

}





