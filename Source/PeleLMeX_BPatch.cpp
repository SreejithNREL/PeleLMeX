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
	ps.get("boundary_dimension", m_boundary_dim);
	if(m_boundary_dim!=0 && m_boundary_dim!=1 && m_boundary_dim!=2)
	{
		amrex::Abort("\nBoundary dimension should be 0,1 or 2");
	}
	ps.get("low_or_high", m_boundary_lo_hi);
	if(m_boundary_lo_hi!=0 && m_boundary_lo_hi!=1)
	{
		amrex::Abort("\nBoundary high low should be 0 or 1");
	}
	if(m_patchtype=="fullboundary")
	{
		//treat it like its a rectangular boundary

	}
	else if(m_patchtype=="circle"){
		ps.get("circle_radius", m_patch_circle_radius);
		if(m_patch_circle_radius<=0.0)
		{
			amrex::Abort("\nPatch type of circle should have radius greater than 0");

		}
		ps.getarr("circle_center", m_patch_circle_center);
	}
	else if(m_patchtype=="rectangle"){
		ps.getarr("patch_rectangle_lo",m_patch_rectangle_lo);
		ps.getarr("patch_rectangle_hi",m_patch_rectangle_hi);

	}
	else if(m_patchtype=="circle-annular"){
		ps.getarr("patch_circ_ann_center",m_patch_circ_ann_center);
		ps.get("patch_circ_ann_inner_radius",m_patch_circ_ann_inner_radius);
		ps.get("patch_circ_ann_outer_radius",m_patch_circ_ann_outer_radius);
		if(m_patch_circ_ann_outer_radius<=0.0 || m_patch_circ_ann_inner_radius<=0.0 || m_patch_circ_ann_inner_radius>=m_patch_circ_ann_outer_radius)
		{
			amrex::Abort("\nRadii do not look good for annular circles");
		}
	}
	else if(m_patchtype=="rectangule-annular"){
		ps.getarr("patch_rect_ann_outer_lo",m_patch_rect_ann_outer_lo);
		ps.getarr("patch_rect_ann_outer_hi",m_patch_rect_ann_outer_hi);
		ps.getarr("patch_rect_ann_inner_lo",m_patch_rect_ann_inner_lo);
		ps.getarr("patch_rect_ann_inner_hi",m_patch_rect_ann_inner_hi);
	}
	else{
		amrex::Abort("\nError! Unknown patch type");
	}

	num_species = ps.countval("species");
	ps.getarr("species", speciesList);

	if (num_species > 0){
		speciesList.resize(num_species);
		speciesIndex.resize(num_species);
		speciesFlux.resize(num_species);
	}
	else
	{
		amrex::Abort("\nError! No species provided to plot flux");
	}

	amrex::Vector<std::string> names;
	pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(names);
	names.resize(names.size());

	for(int n =0;n<names.size();n++){
		speciesIndex[n]=-1;
	}
	for(int m=0;m<num_species;m++){
		for(int n =0;n<names.size();n++){
			if(speciesList[m]==names[n])
			{
				speciesIndex[m]=n;
			}
		}

	}

	for(int n =0;n<num_species;n++){
		if(speciesIndex[n]==-1){
			std::string msg = "\nError! Unable to find species index "+std::to_string(n);
			amrex::Abort(msg);
		}
	}





}


bool BPatch::CheckifPointInside(amrex::GpuArray <amrex::Real, AMREX_SPACEDIM> point_coordinate,const amrex::Geometry& geom) const
{

	bool inside=false;
	const amrex::Real sqrt2 = sqrt(2.0);
	const auto dx = geom.CellSizeArray();
	amrex::Real xp,yp,zp;

#if (AMREX_SPACEDIM == 2)
	xp=point_coordinate[0];
	yp=point_coordinate[1];
#elif (AMREX_SPACEDIM == 3)
	xp=point_coordinate[0];
	yp=point_coordinate[1];
	zp=point_coordinate[2];
#endif


	if(m_patchtype=="circle"){
		if(m_boundary_dim==2)
		{
			amrex::Real rad = sqrt((xp-m_patch_circle_center[0])*(xp-m_patch_circle_center[0])+(yp-m_patch_circle_center[1])*(yp-m_patch_circle_center[1]));
			if(rad<=m_patch_circle_radius+dx[0]*sqrt2/2.0){
				inside=true;
			}

		}

	}
	return inside;

}
