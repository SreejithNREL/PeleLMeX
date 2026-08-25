# AMReX
DIM = 2
COMP = gnu
PRECISION = DOUBLE
USE_EB = FALSE
USE_HYPRE = FALSE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PelePhysics
USE_PLASMA = FALSE
USE_PARTICLES = TRUE
SPRAY_GCM = FALSE
SPRAY_FUEL_NUM = 1

# Optional slurry drying parcel model.
# When TRUE, Exec/Make.PeleLMeX adds -DPELE_SPRAY_DRYING and
# -DPELE_SPRAY_DRYING_NR=$(SPRAY_DRYING_NR)
SPRAY_DRYING = FALSE
SPRAY_DRYING_NR = 20

ifeq ($(SPRAY_GCM), TRUE)
  SpraySuffix = .$(SPRAY_FUEL_NUM)SprayFuel.SprayGCM
else
  SpraySuffix = .$(SPRAY_FUEL_NUM)SprayFuel.SprayMP
endif
ifeq ($(SPRAY_DRYING), TRUE)
  SpraySuffix := $(SpraySuffix).SprayDrying$(SPRAY_DRYING_NR)
endif

USE_MANIFOLD = FALSE
ifeq ($(USE_MANIFOLD), TRUE)
  Eos_Model       = Manifold
  Chemistry_Model = Null
  Transport_Model = Manifold
  Manifold_Dim = $(SPRAY_FUEL_NUM)
  Manifold_Type = Table
  USERSuffix = .Manifold$(Manifold_Dim)d$(Manifold_Type)$(SpraySuffix)
else
  Chemistry_Model = liquid_fuels_nonreacting
  Eos_Model = Fuego
  Transport_Model = Simple
  USERSuffix = .$(Chemistry_Model)$(SpraySuffix)
endif

PELE_HOME ?= ../../..
include $(PELE_HOME)/Exec/Make.PeleLMeX
