include ../make.in

CFLAGS = -O3 -Wall -openmp -w2 -std=gnu99
CFLAGS += -Wunknown-pragmas -Wunused-variable -Wunused-function
CFLAGS += ${OPTFLAGS}

INC = -I./             \
      -I../include \
      -I../external/Yeppp/include

LDFLAGS = -L./

LIBS =  -lifcore

ifeq "${arch}" "mic"
UNIT_TEST_REF_EXE = Test.Ref.MIC
UNIT_TEST_OPT_EXE = Test.Opt.MIC
PERF_TEST_EXE = Test.Perf.MIC
UNIT_TEST_REF_OBJS := mic/testCInt.ref.o mic/screening.ref.o
UNIT_TEST_OPT_OBJS := mic/testCInt.opt.o mic/screening.opt.o
PERF_TEST_OBJS := mic/testCInt2.o mic/screening.o
LIBS_CINT = ../lib/libcint_mic.a \
            ../lib/liboed_mic.a
LIBS_REF = ../lib/liberd_ref_mic.a 
LIBS_OPT = ../lib/liberd_opt_mic.a
else
UNIT_TEST_REF_EXE = Test.Ref
UNIT_TEST_OPT_EXE = Test.Opt
PERF_TEST_EXE = Test.Perf
UNIT_TEST_REF_OBJS := testCInt.ref.o screening.ref.o
UNIT_TEST_OPT_OBJS := testCInt.opt.o screening.opt.o
PERF_TEST_OBJS := testCInt2.o screening.o
LIBS_CINT = ../lib/libcint.a \
            ../lib/liboed.a
LIBS_REF = ../lib/liberd_ref.a
LIBS_OPT = ../lib/liberd_opt.a
endif

all: ${UNIT_TEST_REF_EXE} ${UNIT_TEST_OPT_EXE} ${PERF_TEST_EXE}

${UNIT_TEST_REF_EXE}: $(UNIT_TEST_REF_OBJS) $(LIBS_CINT) ${LIBS_REF}
	$(CC) ${CFLAGS} ${LDFLAGS} $(UNIT_TEST_REF_OBJS) -o $@ ${LIBS} ${LIBS_CINT} ${LIBS_REF}

${UNIT_TEST_OPT_EXE}: $(UNIT_TEST_OPT_OBJS) $(LIBS_CINT) ${LIBS_OPT}
	$(CC) ${CFLAGS} ${LDFLAGS} $(UNIT_TEST_OPT_OBJS) -o $@ ${LIBS} ${LIBS_CINT} ${LIBS_OPT}

${PERF_TEST_EXE}: $(PERF_TEST_OBJS) $(LIBS_CINT) ${LIBS_OPT}
	$(CC) ${CFLAGS} ${LDFLAGS} $(PERF_TEST_OBJS) -o $@ ${LIBS} ${LIBS_CINT} ${LIBS_OPT}

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

%.ref.o : %.c Makefile
	$(CC) ${CFLAGS} -DOPTERD_TEST_REFERENCE ${INC} -c $< -o $@

%.opt.o : %.c Makefile
	$(CC) ${CFLAGS} -DOPTERD_TEST_OPTIMIZED ${INC} -c $< -o $@

mic/%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

mic/%.ref.o : %.c Makefile
	$(CC) ${CFLAGS} -DOPTERD_TEST_REFERENCE ${INC} -c $< -o $@

mic/%.opt.o : %.c Makefile
	$(CC) ${CFLAGS} -DOPTERD_TEST_OPTIMIZED ${INC} -c $< -o $@

clean:
	rm -f *.o *.s *.d *~ mic/*.o testCInt*.o Test.*
