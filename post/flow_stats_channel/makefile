
#F90	= $(FC)

FC	= gfortran
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3
FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-200 -fno-align-commons -fbounds-check -cpp  #-std=f2008ts -g
#FOPTS	= -align none -mcmodel medium -warn all
FFLAGS	= $(FOPTS) $(OLEVEL)

LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=

FLWOBJS = \
./module.f90 \
./main.f90 \
./read_fields.f90 

MOD	=
OBJS	= $(FLWOBJS)
EXEC    =  ./stats

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) -o $@

clean:
	$(RM) $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
