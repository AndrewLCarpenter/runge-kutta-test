.SUFFIXES: .o .c .f90 .F90 .mod

CC = gcc
#FC = ifort
 FC = gfortran

#FCFLAGS = -Wuninitialized -Wunused -pedantic -Wall -Wsurprising -fbacktrace -fbounds-check -W -O1 $(INCLUDESDIR)
#FCFLAGS = -Wuninitialized -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
#FCFLAGS = -check uninit -warn all -no-wrap-margin -WB -O3 $(INCLUDESDIR)
#FCFLAGS = -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
 FCFLAGS = -Wuninitialized -Wsurprising -O3 $(INCLUDESDIR)
  CFLAGS = -Ofast  $(INCLUDESDIR)

SRCS = precision_vars.f90\
       unary_mod.f90\
       eispack_module.f90\
       blas_module.f90\
       OTD_Module.f90\
       SBP_Coef_Module.f90\
       Jacobian_CSR_Mod.f90\
       control_variables.f90\
       poly_fit_Mod.f90\
       runge_kutta.f90\
       Fourier_Coef_Module.f90\
       Stage_Value_module.f90\
       QR_Module.f90\
       matvec_module.f90\
       ilut_module.f90\
       Kuramoto_Sivashinsky.f90\
       Boldrighini_Franceschini.f90\
       Burgers_Mod.f90\
       vanderPol.f90\
       Pareschi.f90\
       Kaps.f90\
       Kreiss.f90\
       Lorenz.f90\
       Rossler_Chaos.f90\
       Oregonator.f90\
       brusselator_2.f90\
       Boscarino31_mod.f90\
       Broadwell_mod.f90\
       Charney_DeVore6.f90\
       Charney_DeVore10.f90\
       problemsub.f90\
       output_module.f90\
       Newton.f90\
       time_integ.f90\
       test_cases.f90


      
#OBJS2 = $(SRCS:.F90=.o)

#OBJS= $(OBJS2) \
$(OBJS1)

OBJS = precision_vars.o\
       unary_mod.o\
       eispack_module.o\
       blas_module.o\
       OTD_Module.o\
       SBP_Coef_Module.o\
       Jacobian_CSR_Mod.o\
       control_variables.o\
       poly_fit_Mod.o\
       runge_kutta.o\
       Fourier_Coef_Module.o\
       Stage_Value_module.o\
       QR_Module.o\
       matvec_module.o\
       ilut_module.o\
       Kuramoto_Sivashinsky.o\
       Boldrighini_Franceschini.o\
       Burgers_Mod.o\
       vanderPol.o\
       Pareschi.o\
       Kaps.o\
       Kreiss.o\
       Lorenz.o\
       Rossler_Chaos.o\
       Oregonator.o\
       brusselator_2.o\
       Boscarino31_mod.o\
       Broadwell_mod.o\
       Charney_DeVore6.o\
       Charney_DeVore10.o\
       problemsub.o\
       output_module.o\
       Newton.o\
       time_integ.o\
       test_cases.o

      
#LINK = gfortran

TARGET = test_cases90.out

.PHONY: depend clean

all: $(TARGET)
	@echo  Runge-Kutta compiled. Executable in test_cases90.out

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
