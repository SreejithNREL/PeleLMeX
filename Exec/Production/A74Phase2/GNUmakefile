TOP = ../../..
PELELMEX_HOME      ?= ${TOP}
AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex
PELE_PHYSICS_HOME  ?= ${PELELMEX_HOME}/Submodules/PelePhysics
AMREX_HYDRO_HOME   ?= ${PELELMEX_HOME}/Submodules/AMReX-Hydro
PELEMP_HOME        ?= ${PELELMEX_HOME}/Submodules/PeleMP
#HYPRE_HOME         = /ccs/home/lesclapez/hypre/install-2807


# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE
USE_EB          = FALSE
USE_HYPRE       = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

#PELE_USE_MAGMA  = TRUE

# PeleLMeX
CEXE_sources   += PeleLMeX_DeriveUserDefined.cpp

# PeleMP
USE_PARTICLES   = TRUE
SPRAY_FUEL_NUM  = 1
CEXE_sources   += SprayParticlesInitInsert.cpp
PELELM_USE_SPRAY = TRUE

# PelePhysics
Chemistry_Model = dodmethair_4sp
Eos_Model	= Fuego
Transport_Model = Simple

include $(PELELMEX_HOME)/Exec/Make.PeleLMeX
