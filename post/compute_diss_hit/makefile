#F90	= $(FC)

FC	= gfortran
LD	= $(FC)
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3 
FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp  #-std=f2008ts -g 
FFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=

FLWOBJS = \
./module.f90 \
./main.f90 \
./read_fields.f90 


MOD	= 
OBJS	= $(FLWOBJS) 
EXEC    =  ./epsilon

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS)  -o $@ 

clean:
	$(RM) $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
