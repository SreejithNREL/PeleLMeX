#ifndef BPATCH_H
#define BPATCH_H

#include <PeleLMeX.H>

class BPatch
{
public:

	BPatch() = default;
	BPatch(const std::string& patch_name, const amrex::Geometry& geom);
	bool CheckifPointInside(amrex::GpuArray <amrex::Real, AMREX_SPACEDIM> ,const amrex::Geometry& ) const;
	int getPatchBoundaryDim(){ return m_boundary_dim;}
	int getPatchBoundaryLoHi(){return m_boundary_lo_hi;}
	int getNumSpecies(){return num_species;}
	int getSpeciesIdx(int i){return speciesIndex[i];}
	void setSpeciesFlux(int i,amrex::Real flux){speciesFlux[i]=flux;}



protected:
	std::string m_patchname;
	std::string m_patchtype;	//values: "fullboundary","circle","rectangle","annular-circle","annular-rectangle"
	int m_boundary_dim;
	int m_boundary_lo_hi;
	int num_species;
	amrex::Vector<std::string> speciesList;
	amrex::Vector<int> speciesIndex;
	amrex::Vector<amrex::Real> speciesFlux;

	//circle variables
	amrex::Vector<amrex::Real> m_patch_circle_center;
	amrex::Real  m_patch_circle_radius;

	//rectangle variables
	amrex::Vector<amrex::Real> m_patch_rectangle_lo;
	amrex::Vector<amrex::Real> m_patch_rectangle_hi;

	//circle annular variables
	amrex::Vector<amrex::Real> m_patch_circ_ann_center;
	amrex::Real  m_patch_circ_ann_inner_radius;
	amrex::Real  m_patch_circ_ann_outer_radius;

	//Rectangular annular variables
	amrex::Vector<amrex::Real> m_patch_rect_ann_outer_lo;
	amrex::Vector<amrex::Real> m_patch_rect_ann_outer_hi;
	amrex::Vector<amrex::Real> m_patch_rect_ann_inner_lo;
	amrex::Vector<amrex::Real> m_patch_rect_ann_inner_hi;

};

#endif
