#ifndef BPATCH_H
#define BPATCH_H

#include <PeleLMeX.H>

class BPatch
{

public:

	BPatch()=default;
	~BPatch() = default;

	BPatch(const std::string& patch_name, const amrex::Geometry& geom);

	struct BpatchDataContainer
	{
		int m_boundary_dim;
		int m_boundary_lo_hi;
		int num_species;
		int m_patchtype_num;
		int * speciesIndex=nullptr;
		amrex::Real * speciesFlux=nullptr;

		#if (AMREX_SPACEDIM == 1)
			//declare nothing
		#elif (AMREX_SPACEDIM == 2)
			//In 2D, a boundary patch is a line
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_line_center;
			amrex::Real  m_patch_line_radius;
		#else
			//circle variables
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_circle_center;
			amrex::Real  m_patch_circle_radius;

			//rectangle variables
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_rectangle_lo;
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_rectangle_hi;

			//circle annular variables
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_circ_ann_center;
			amrex::Real  m_patch_circ_ann_inner_radius;
			amrex::Real  m_patch_circ_ann_outer_radius;

			//Rectangular annular variables
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_rect_ann_outer_lo;
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_rect_ann_outer_hi;
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_rect_ann_inner_lo;
			amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_patch_rect_ann_inner_hi;
		#endif

		AMREX_GPU_DEVICE
		AMREX_FORCE_INLINE
		bool CheckifPointInside(amrex::GpuArray <amrex::Real, AMREX_SPACEDIM> point_coordinate,
				amrex::Real dx) const
		{
			bool inside=false;
			const amrex::Real sqrt2 = sqrt(2.0);
			//auto dx = geom.CellSizeArray();
			amrex::Real xp,yp,zp;

			#if (AMREX_SPACEDIM == 2)
				xp=point_coordinate[0];
				yp=point_coordinate[1];
			#elif (AMREX_SPACEDIM == 3)
				xp=point_coordinate[0];
				yp=point_coordinate[1];
				zp=point_coordinate[2];
			#endif

			//if(m_patchtype=="circle"){
				if(m_boundary_dim==2)
				{
					amrex::Real rad = sqrt((xp-m_patch_circle_center[0])*(xp-m_patch_circle_center[0])+(yp-m_patch_circle_center[1])*(yp-m_patch_circle_center[1]));
					if(rad<=m_patch_circle_radius+dx*sqrt2/2.0){
						inside=true;
					}
				}
			//}
			return inside;
		}

	};

	void allocate()
	 {
	   if (!m_device_allocated) {
		   m_bpdata_d = (BpatchDataContainer*)amrex::The_Arena()->alloc(sizeof(m_bpdata_h));
	     m_device_allocated = true;
	     sync_to_device();
	   }
	 }

	void deallocate()
	{
	   if (m_host_allocated) {
	     amrex::The_Pinned_Arena()->free(m_bpdata_h.speciesIndex);
	     amrex::The_Pinned_Arena()->free(m_bpdata_h.speciesFlux);
	   }
	   if (m_device_allocated) {
	     amrex::The_Arena()->free(m_bpdata_d);
	   }
	 }

	  void sync_to_device()
	  {
	    if (!m_device_allocated) {
	      amrex::Abort("Device params not allocated yet");
	    } else {
	      amrex::Gpu::copy(
	        amrex::Gpu::hostToDevice, &m_bpdata_h, &m_bpdata_h + 1, m_bpdata_d);
	    }
	  }

	// Accessors
	BpatchDataContainer getHostData() { return m_bpdata_h; }
	BpatchDataContainer* getDeviceData() { return m_bpdata_d; }

	std::string m_patchname;
	std::string m_patchtype;	//values: "fullboundary","circle","rectangle","annular-circle","annular-rectangle"
	amrex::Vector<std::string> speciesList;

protected:
	BpatchDataContainer m_bpdata_h;
	BpatchDataContainer* m_bpdata_d = nullptr;
	bool m_device_allocated{false};
	bool m_host_allocated{false};
};

#endif
