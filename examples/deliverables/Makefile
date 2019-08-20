#BHEADER**********************************************************************
#
# Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. Written by 
# Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
# Dobrev, et al. LLNL-CODE-660355. All rights reserved.
# 
# This file is part of XBraid. For support, post issues to the XBraid Github page.
# 
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2.1 dated February 1999.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with this program; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#EHEADER**********************************************************************

##################################################################
# Import machine specific compilers, options, flags, etc.. 
##################################################################

BRAID_DIR=../../braid
include ../../makefile.inc

##################################################################
# Build exmaples 
##################################################################

HYPRE_DIR = ../../../hypre/src/hypre
HYPRE_FLAGS = -I$(HYPRE_DIR)/include
HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE
HYPRE_LIB_FILE = $(HYPRE_DIR)/lib/libHYPRE.a

MFEM_DIR = ../../../mfem
MFEM_CONFIG_MK = $(MFEM_DIR)/config/config.mk
MFEM_LIB_FILE = mfem_is_not_built
-include $(MFEM_CONFIG_MK)

BRAID_FLAGS = -I$(BRAID_DIR)
BRAID_LIB_FILE = $(BRAID_DIR)/libbraid.a

C_NOHYPRE = ex-01 ex-01-adjoint ex-01-optimization ex-01-refinement ex-01-expanded ex-01-expanded-bdf2 ex-02 ex-04 ex-04-serial ex-04-omgrit
CPP_NOHYPRE = ex-01-pp 
F_NOHYPRE = ex-01-expanded-f
C_EXAMPLES = ex-03 ex-03-serial
# Note: .cpp examples will be linked with mfem
#CXX_EXAMPLES = ex-04

.PHONY: all clean cleanout

.SUFFIXES:
.SUFFIXES: .c .cpp

# put this rule first so it becomes the default
all: $(C_NOHYPRE) $(CPP_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES)


advec-diff-implicit: advec-diff-implicit.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

advec-diff-rms: advec-diff-rms.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

advec-diff-upwind: advec-diff-upwind.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

advec-diff-serial: advec-diff-serial.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

visc-burgers-explicit: visc-burgers-explicit.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

visc-burgers-implicit: visc-burgers-implicit.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

visc-burgers-serial: visc-burgers-serial.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)


# Rule for compiling .c files
%: %.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(HYPRE_FLAGS) $(@).c -o $@\
 $(BRAID_LIB_FILE) $(HYPRE_LIB) $(LFLAGS)

# Rule for compiling .cpp files; links with mfem
%: %.cpp $(BRAID_LIB_FILE) $(MFEM_LIB_FILE) $(MFEM_CONFIG_MK)
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) $(BRAID_FLAGS) $(MFEM_FLAGS) \
	$< -o $@ $(MFEM_LIBS) $(BRAID_LIB_FILE) $(LFLAGS)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

# Generate an error message if the Hypre library is not built and exit
$(HYPRE_LIB_FILE):
	$(error The Hypre library is not built, unable to build ex-03)

clean: cleanout
	rm -f *.o $(C_NOHYPRE) $(CPP_NOHYPRE) $(F_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES) $(F_EXAMPLES) *ror_norm* *_err_* *_mesh* *_sol_*
	rm -rf *.dSYM

cleanout:
	rm -f ex*.out.*

