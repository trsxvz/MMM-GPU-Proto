# ROCm from its Python packages (TheRock), integrated GPU gfx1036
ROCM_VENV ?= $(HOME)/Bureau/rocm-sdk-therock/venv
ROCM ?= $(shell $(ROCM_VENV)/bin/rocm-sdk path --root)
ARCH ?= gfx1036
# MFEM compiled with HIP for gfx1036 (scripts/build-hip-amd.sh)
MFEM_INSTALL ?= ../../deps/hip-amd/install/mfem-gfx1036-hypre-hip

all: idle-loss

# plain HIP
idle-loss: idle-loss.cpp
	$(ROCM)/bin/hipcc --offload-arch=$(ARCH) -isystem $(ROCM)/include idle-loss.cpp -o idle-loss

# through MFEM, as first seen in mfem-mgis
lost-kernel: lost-kernel.cpp
	$(MAKE) -f Makefile.mfem MFEM_INSTALL=$(MFEM_INSTALL)

.PHONY: all lost-kernel
