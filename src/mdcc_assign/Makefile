EXECUTABLES = mdcc_assign
CC=g++
ICCFLAGS = -O3 #-w2
LAPACKLIB = -lmkl_lapack95_ilp64
PATH_LAPACKLIB = /opt/intel/mkl/lib/intel64/
#LAPACKLIB = -llapack
#PATH_LAPACKLIB = ${HOME}/lib
BOOSTINC = $(HOME)/include

OBJS = PliGauss.o Config.o Read.o Write.o\
       Inact.o Type.o Coord.o \
       Gaussian.o math.o \
       main.o

#OBJS =$(SRCS:.c=.o)

.SUFFIXES: .o .cpp

$(EXECUTABLES) : $(OBJS)
	rm -f $(EXECUTABLES)
	$(CC) -o $(EXECUTABLES) $(OBJS) -L$(PATH_LAPACKLIB) ${LAPACKLIB} --static
.cpp.o:
	$(CC) -c $(ICCFLAGS) -I$(BOOSTINC)  $<

$(OBJS): define.h
clean :
	-rm -f *~ *.o $(EXECUTABLES)
