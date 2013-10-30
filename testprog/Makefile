include ../make.in

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
UNIT_TEST_REF_EXE = Test.Ref.MIC
UNIT_TEST_OPT_EXE = Test.Opt.MIC
PERF_TEST_EXE = Test.Perf.MIC
UNIT_TEST_REF_OBJS := mic/testCInt.ref.o
UNIT_TEST_OPT_OBJS := mic/testCInt.opt.o
PERF_TEST_OBJS := mic/testCInt2.o
LIBS_CINT = ../lib/libcint_mic.a \
            ../lib/liberd_mic.a  \
            ../lib/liboed_mic.a
else
UNIT_TEST_REF_EXE = Test.Ref
UNIT_TEST_OPT_EXE = Test.Opt
PERF_TEST_EXE = Test.Perf
UNIT_TEST_REF_OBJS := testCInt.ref.o
UNIT_TEST_OPT_OBJS := testCInt.opt.o
PERF_TEST_OBJS := testCInt2.o
LIBS_CINT = ../lib/libcint.a \
            ../lib/liboed.a
endif

all: ${UNIT_TEST_REF_EXE} ${UNIT_TEST_OPT_EXE} ${PERF_TEST_EXE}

${UNIT_TEST_REF_EXE}: $(UNIT_TEST_REF_OBJS) $(LIBS_CINT) ../lib/liberd_ref.a
	$(CC) ${CFLAGS} ${LDFLAGS} $(UNIT_TEST_REF_OBJS) -o $@ ${LIBS} ${LIBS_CINT} ../lib/liberd_ref.a

${UNIT_TEST_OPT_EXE}: $(UNIT_TEST_OPT_OBJS) $(LIBS_CINT) ../lib/liberd_opt.a
	$(CC) ${CFLAGS} ${LDFLAGS} $(UNIT_TEST_OPT_OBJS) -o $@ ${LIBS} ${LIBS_CINT} ../lib/liberd_opt.a

${PERF_TEST_EXE}: $(PERF_TEST_OBJS) $(LIBS_CINT) ../lib/liberd_opt.a
	$(CC) ${CFLAGS} ${LDFLAGS} $(PERF_TEST_OBJS) -o $@ ${LIBS} ${LIBS_CINT} ../lib/liberd_opt.a

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

%.ref.o : %.c Makefile
	$(CC) ${CFLAGS} -DOPTERD_TEST_REFERENCE ${INC} -c $< -o $@

%.opt.o : %.c Makefile
	$(CC) ${CFLAGS} -DOPTERD_TEST_OPTIMIZED ${INC} -c $< -o $@

mic/%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

clean:
	rm -f *.o *.s *.d *~ mic/* testCInt*.o testPFock
