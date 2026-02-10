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
get_interface.f90 \
main.f90

MOD	= 



OBJS	= $(FLWOBJS)
EXEC    =  ./intearea

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS)  -o $@

clean:
	$(RM) $(EXEC) *.mod *.o

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
