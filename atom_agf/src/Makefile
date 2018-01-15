# Makefile for EPW

include make.sys

#
# use recursive definitions in Makefile to modify the PATH for modules
#

#IFLAGS   = -I../../include  
#MODFLAGS = -I./ -I../../Modules -I../../iotk/src \
#           -I../../PW -I../../PH -I../../PP
#LIBOBJS  = ../../flib/ptools.a ../../flib/flib.a \
#           ../../clib/clib.a ../../iotk/src/libiotk.a

OBJS = \
avgmods.o \
avgcalc.o 
#histcalc.o


default : green

all : green

green : $(OBJS)
	$(LD) -o green.x \
	$(OBJS) $(LDFLAGS) $(LIBS)
	- ( mkdir -p ../bin ; cd ../bin ; ln -fs ../src/green.x . )

clean :
	- /bin/rm -f  *.o *~ *.F90 *.d *.mod *.i 

