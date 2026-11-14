# Built against the installation of MFEM of the prototype (MFEM_DIR).
MFEM_DIR ?= /lus/work/RES1/genden15/tchenaille/mmm/install/mfem
include $(MFEM_DIR)/share/mfem/config.mk
partition-mesh: partition_mesh.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)
