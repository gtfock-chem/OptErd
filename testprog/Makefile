include ../make.in

OBJS_PFOCK := testPFock.o
OBJS_SCR := testScr.o

CFLAGS = -O3 -Wall -openmp -w2
CFLAGS += -Wunknown-pragmas -Wunused-variable #-Wunused-function
CFLAGS += ${OPTFLAGS}

INC = -I./             \
      -I../include     \
      -I${MPI_INCDIR}  \
      -I${BLAS_INCDIR} \
      -I${GA_INCDIR}

LDFLAGS = -L./ -L${MPI_LIBDIR} -L${BLAS_LIBDIR}

LIBS =  -lifcore

ifeq "${arch}" "mic"
EXE1 = testCInt_mic
EXE2 = testCInt2_mic
OBJS_CINT := mic/testCInt.o
OBJS_CINT2 := mic/testCInt2.o
LIBS_CINT = ../lib/libcint_mic.a \
            ../lib/liberd_mic.a  \
            ../lib/liboed_mic.a
else
EXE1 = testCInt
EXE2 = testCInt2
OBJS_CINT := testCInt.o
OBJS_CINT2 := testCInt2.o
LIBS_CINT = ../lib/libcint.a \
            ../lib/liberd.a  \
            ../lib/liboed.a
endif

LIBS_PFOCK = ../lib/libpfock.a ${MPI_LIBS} ${SCALAPACK_LIBS} ${GA_LIBS} ${BLAS_LIBS}

all: ${EXE1} ${EXE2}

${EXE1}: $(OBJS_CINT)
	$(CC) ${CFLAGS} ${LDFLAGS} $(OBJS_CINT) -o $@ ${LIBS} ${LIBS_CINT}

${EXE2}: $(OBJS_CINT2)
	$(CC) ${CFLAGS} ${LDFLAGS} $(OBJS_CINT2) -o $@ ${LIBS} ${LIBS_CINT}

testPFock: $(OBJS_PFOCK)
	$(MPICC) ${CFLAGS} ${LDFLAGS} $(OBJS_PFOCK) -o $@ ${LIBS} ${LIBS_PFOCK} ${LIBS_CINT}

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@ 

mic/%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

clean:
	rm -f *.o *.s *.d *~ mic/* testCInt*.o testPFock
