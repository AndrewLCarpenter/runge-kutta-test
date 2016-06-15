.SUFFIXES: .o .c .f90 .F90

CC = gcc
FC = gfortran

FCFLAGS = -Wmaybe-uninitialized -Wsurprising -fbacktrace -fbounds-check -O3 $(INCLUDESDIR)
CFLAGS = -Ofast  $(INCLUDESDIR)

SRCS = SBP_Coef_Module.f90\
       driver.f90

OBJS = SBP_Coef_Module.o\
       driver.o

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

