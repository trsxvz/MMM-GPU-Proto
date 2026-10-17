# builds the example ex1p of PETSc of MFEM 4.10 (unmodified copy of
# examples/petsc/ex1p.cpp) against an installation of MFEM built with CUDA
# and PETSc: make MFEM_DIR=<installation of MFEM>
ifndef MFEM_DIR
$(error give MFEM_DIR=<installation of MFEM>)
endif
include $(MFEM_DIR)/share/mfem/config.mk

ex1p: ex1p.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

clean:
	rm -f ex1p
