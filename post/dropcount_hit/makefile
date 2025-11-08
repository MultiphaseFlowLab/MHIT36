FC	= gfortran
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
# OLEVEL	= -O2
OLEVEL	= -g -fbacktrace  -fbounds-check -Waliasing -Wunderflow -Wsurprising -fbacktrace -fcheck=all -Wall -fcheck=all
FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp  #-std=f2008ts -g

FFLAGS	= $(FOPTS) $(OLEVEL)

#LIBS    =  -I/usr/include/ -I/usr/local/include -L$FFTW_LIB -lfftw3 -lm -ldl


LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=

FLWOBJS = \
module.f90 \
read_input.f90 \
read_phi.f90 \
calculate_deq.f90 \
flood_fill.f90 \
get_interface.f90 \
main.f90
# ./read_fields.f90 \
# ./spectral_to_phys.f90 \
# ./dctz_bwd.f90 \
# ./fftx_bwd.f90 \
# ./ffty_bwd.f90 \
# ./create_plan.f90 \
# ./destroy_plan.f90 \
# ./get_interface.f90\
# ./flood_fill.f90 

MOD	= 



OBJS	= $(FLWOBJS)
EXEC    =  ./DropCountSize.x

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS)  -o $@

clean:
	$(RM) $(EXEC) *.mod *.o

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
