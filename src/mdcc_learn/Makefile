EXECUTABLES = mdcc_learn
#CC=g++
CC=icpc
##

LAPACKLIB = -lmkl_lapack95_ilp64
# PATH_LAPACKLIB = /opt/intel/mkl/lib/intel64/
# LAPACKLIB = -llapack
# PATH_LAPACKLIB = ${HOME}/lib
ICCFLAGS = -ansi -O3 #-g #-w1 

OBJS = Kprml.o Config.o Read.o Write.o\
       Coord.o\
       FeatureDefinition.o \
       EMAlgorithm.o \
       VBgmm.o \
       VBfull.o \
       Gaussian.o \
       math.o \
       main.o

#OBJS =$(SRCS:.c=.o)

.SUFFIXES: .o .cpp

$(EXECUTABLES) : $(OBJS)
	rm -f $(EXECUTABLES)
	$(CC) -o $(EXECUTABLES) $(OBJS) -L$(PATH_LAPACKLIB) ${LAPACKLIB} # --static
.cpp.o:
	$(CC) -c $(ICCFLAGS)  $<

$(OBJS): define.h
clean :
	-rm -f *~ *.o $(EXECUTABLES)
