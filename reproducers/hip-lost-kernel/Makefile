# MFEM compiled with HIP for gfx1036 (scripts/build-hip-amd.sh)
MFEM_INSTALL ?= ../../deps/hip-amd/install/mfem-gfx1036-hypre-hip
include $(MFEM_INSTALL)/share/mfem/config.mk
lost-kernel: lost-kernel.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) lost-kernel.cpp -o lost-kernel $(MFEM_LIBS)
