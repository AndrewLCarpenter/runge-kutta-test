.SUFFIXES: .o .c .f90 .F90 .mod

CC = gcc
FC = gfortran

#FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
FCFLAGS = -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)

SRCS = test_cases.f90\
       JACOB2.f90\
       vanderPol.f90\
       Pureschi.f90\
       Kaps.f90\
       Kreiss.f90\
       problemsub.f90\
       fit.f90\
       gser.f90\
       gammln.f90\
       gammq.f90\
       gcf.f90\
       rungeadd.f90
       #Allocate_CSR_Storage.f90\
       #CSR_variables.f90\
       #Invert_Jacobian.f90\
       #Jacobian.f90\
       #Jacobian_CSR.f90\
       #precision_vars.f90


#OBJS2 = $(SRCS:.F90=.o)

#OBJS= $(OBJS2) \
$(OBJS1)

OBJS = test_cases.o\
       JACOB2.o\
       vanderPol.o\
       Pureschi.o\
       Kaps.o\
       Kreiss.o\
       problemsub.o\
       fit.o\
       gser.o\
       gammln.o\
       gammq.o\
       gcf.o\
       rungeadd.o
       #Allocate_CSR_Storage.o\
       #CSR_variables.o\
       #Invrt_Jacobian.o\
       #Jacobian.o\
       #precision_vars.o\
       #Jacobian_CSR.o
      
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
