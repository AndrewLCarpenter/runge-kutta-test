.SUFFIXES: .o .c .f90 .F90 .mod

CC = gcc
FC = ifort
#FC = gfortran

#FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
# FCFLAGS = -Wuninitialized -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
FCFLAGS = -check uninit -WB -O3 $(INCLUDESDIR)
#FCFLAGS = -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)

SRCS = precision_vars.f90\
       SBP_Coef_Module.f90\
       control_variables.f90\
       poly_fit_Mod.f90\
       runge_kutta.f90\
       output_module.f90\
       Stage_Value_module.f90\
       QR_Module.f90\
       Jacobian_CSR_Mod.f90\
       matvec_module.f90\
       blas_module.f90\
       eispack_module.f90\
       ilut_module.f90\
       unary_mod.f90\
       Burgers_Mod.f90\
       vanderPol.f90\
       Pureschi.f90\
       Kaps.f90\
       Kreiss.f90\
       Lorenz.f90\
       Rossler_Chaos.f90\
       Oregonator.f90\
       brusselator.f90\
       problemsub.f90\
       Newton.f90\
       test_cases.f90


      
#OBJS2 = $(SRCS:.F90=.o)

#OBJS= $(OBJS2) \
$(OBJS1)

OBJS = precision_vars.o\
       SBP_Coef_Module.o\
       control_variables.o\
       poly_fit_Mod.o\
       runge_kutta.o\
       output_module.o\
       Stage_Value_module.o\
       QR_Module.o\
       Jacobian_CSR_Mod.o\
       matvec_module.o\
       blas_module.o\
       eispack_module.o\
       ilut_module.o\
       unary_mod.o\
       Burgers_Mod.o\
       problemsub.o\
       vanderPol.o\
       Pureschi.o\
       Kaps.o\
       Kreiss.o\
       Lorenz.o\
       Rossler_Chaos.o\
       Oregonator.o\
       brusselator.o\
       Newton.o\
       test_cases.o

      
#LINK = gfortran

TARGET = test_cases90.out

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
