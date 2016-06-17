.SUFFIXES: .o .c .f90 .F90

CC = gcc
FC = gfortran

FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)

SRCS = precision_vars.f90\
       blas_module.f90\
       matvec_module.f90\
       eispack_module.f90\
       ilut_module.f90\
       unary_mod.f90\
       SBP_Coef_Module.f90\
       Driver.f90

OBJS = precision_vars.o\
       blas_module.o\
       matvec_module.o\
       eispack_module.o\
       ilut_module.o\
       unary_mod.o\
       SBP_Coef_Module.o\
       Driver.o

TARGET = driver

.PHONY: depend clean

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

