# Compiler
FC = mpif90
LD = $(FC)
 
# Paths (Modify if necessary)
ROOT_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
$(info $$ROOT_DIR is [${ROOT_DIR}])
PARENT_DIR := $(abspath $(ROOT_DIR)/..)
$(info $$PARENT_DIR is [${PARENT_DIR}])
CUDECOMP_DIR = $(PARENT_DIR)/cuDecomp/build
CUDA_DIR = /leonardo/prod/opt/compilers/cuda/12.3/none
EXTRA_DIR = /leonardo/prod/opt/compilers/cuda/12.3/none/compat

export $(CUDECOMP_DIR)/lib

# Compiler and Linker Flags
FFLAGS = -fast -acc -gpu=mem:managed -Minfo=accel -Mfree -Mpreprocess -cpp -cuda  -I$(CUDECOMP_DIR)/include/
LIBS = -L$(CUDECOMP_DIR)/lib/ -L$(CUDA_DIR)/lib64/ -lcudecomp_fort -lcudecomp -cudalib=cufft  -lnvToolsExt

# Source and Object Files
OBJS = module.f90 main.f90 readinput.f90 readwrite.f90
EXEC = mhit36

$(EXEC):        $(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
