.SUFFIXES: .o .c .f90 .F90

CC = gcc
FC = mpif90

FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)

SRCS = precision_vars.f90\
       datatypes.f90\
			 aerodynamic_coefficients.f90\
			 error_bc_no_slip_wall.f90\
			 error_heat_entropy_flow_wall_bc.f90\
			 timeinteg.f90\
			 polyinit.f90\
			 physicsindependent.f90\
			 ssdcsolver.f90 

#OBJS2 = $(SRCS:.F90=.o)

#OBJS= $(OBJS2) \
$(OBJS1)

OBJS = precision_vars.o\
			 datatypes.o\
			 variables.o\
			 referencevariables.o\
			 controlvariables.o\
			 error_bc_no_slip_wall.o\
			 error_heat_entropy_flow_wall_bc.o\
			 timeinteg.o\
			 polyinit.o\
			 physicsindependent.o\
			 ssdcsolver.o\
			 metiselementpartition.o

#LINK = gfortran

TARGET = SSDCNSE

.PHONY: depend clean

#metiselementpartition.o: metiselementpartition.c


all: $(TARGET)
	@echo  Simple compiler named SSDCNSE has been compiled


$(TARGET) : $(OBJS)
	$(FC) $(FCFLAGS) $(INCLUDESDIR) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) -c $< -o $@

%.o : %.F90
	$(FC) -c $(FCFLAGS) -c $< -o $@

%.o : %.c
	$(CC) -c $(CFLAGS) -c $< -o $@


clean:
	rm -f *.o *.mod *~ $(TARGET)

depend: $(SRCS)
	makedepend $(INCLUDESDIR) $^

# DO NOT DELETE THIS LINE -- make depend needs it				
