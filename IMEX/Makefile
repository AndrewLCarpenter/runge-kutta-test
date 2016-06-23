.SUFFIXES: .o .c .f90 .F90 .mod

CC = gcc
FC = gfortran

#FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
FCFLAGS = -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)

SRCS = precision_vars.f90\
       poly_fit_Mod.f90\
       CSR_Variables.f90\
       output_module.f90\
       Invert_Jacobian.f90\
       Allocate_CSR_Storage.f90\
       JACOB2.f90\
       Jacobian.f90\
       Newton_Iteration.f90\
       Stage_Value_Predictor.f90\
       rungeadd.f90\
       test_cases.f90\
       vanderPol.f90\
       Pureschi.f90\
       Kaps.f90\
       Kreiss.f90\
       Lorenz.f90\
       problemsub.f90
      
       


#OBJS2 = $(SRCS:.F90=.o)

#OBJS= $(OBJS2) \
$(OBJS1)

OBJS = precision_vars.o\
       poly_fit_Mod.o\
       CSR_Variables.o\
       output_module.o\
       Invert_Jacobian.o\
       Allocate_CSR_Storage.o\
       Jacobian_CSR.o\
       JACOB2.o\
       Jacobian.o\
       Newton_Iteration.o\
       Stage_Value_Predictor.o\
       rungeadd.o\
       test_cases.o\
       vanderPol.o\
       Pureschi.o\
       Kaps.o\
       Kreiss.o\
       Lorenz.o\
       problemsub.o
      
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
