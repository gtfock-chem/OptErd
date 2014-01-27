include ../make.in

CFLAGS = -O3 -Wall -openmp -w2 -std=gnu99
CFLAGS += -Wunknown-pragmas -Wunused-variable -Wunused-function
CFLAGS += ${OPTFLAGS}

INC = -I./  \
      -I../external/Yeppp/include

LDFLAGS = -L./

LIBS =  -lifcore

ifeq "${arch}" "mic"
UNIT_TEST_REF_EXE = Test.Ref.MIC
UNIT_TEST_OPT_EXE = Test.Opt.MIC
PERF_TEST_EXE = Test.Perf.MIC
LEGACY_PERF_TEST_EXE = Test.Perf.legacy.MIC
UNIT_TEST_REF_OBJS := mic/testCInt.ref.o mic/screening.ref.o
UNIT_TEST_OPT_OBJS := mic/testCInt.opt.o mic/screening.opt.o
PERF_TEST_OBJS := mic/testCInt2.opt.o mic/screening.opt.o
LEGACY_PERF_TEST_OBJS := mic/testCInt2.ref.o mic/screening.ref.o
LIBS_CINT_REF = ../legacy/lib/libcint_mic.a \
                ../legacy/lib/liboed_mic.a
LIBS_CINT_OPT = ../lib/libcint_mic.a \
                ../lib/liboed_mic.a
LIBS_REF = ${LIBS_CINT_REF}
LIBS_OPT = ${LIBS_CINT_OPT}
LIBS_REF += ../legacy/lib/liberd_mic.a
LIBS_OPT += ../lib/liberd_mic.a

else
UNIT_TEST_REF_EXE = Test.Ref
UNIT_TEST_OPT_EXE = Test.Opt
PERF_TEST_EXE = Test.Perf
LEGACY_PERF_TEST_EXE = Test.Perf.legacy
UNIT_TEST_REF_OBJS := testCInt.ref.o screening.ref.o
UNIT_TEST_OPT_OBJS := testCInt.opt.o screening.opt.o
PERF_TEST_OBJS := testCInt2.opt.o screening.opt.o
LEGACY_PERF_TEST_OBJS := testCInt2.ref.o screening.ref.o
LIBS_CINT_REF = ../legacy/lib/libcint.a \
                ../legacy/lib/liboed.a
LIBS_CINT_OPT = ../lib/libcint.a \
                ../lib/liboed.a
LIBS_REF = ${LIBS_CINT_REF}
LIBS_OPT = ${LIBS_CINT_OPT}
LIBS_REF += ../legacy/lib/liberd.a
LIBS_OPT += ../lib/liberd.a
endif


all: ${UNIT_TEST_REF_EXE} ${UNIT_TEST_OPT_EXE} ${PERF_TEST_EXE} ${LEGACY_PERF_TEST_EXE}

${UNIT_TEST_REF_EXE}: $(UNIT_TEST_REF_OBJS) ${LIBS_REF}
	$(CC) ${CFLAGS} ${LDFLAGS} $(UNIT_TEST_REF_OBJS) -o $@ ${LIBS} ${LIBS_REF}

${UNIT_TEST_OPT_EXE}: $(UNIT_TEST_OPT_OBJS) ${LIBS_OPT}
	$(CC) ${CFLAGS} ${LDFLAGS} $(UNIT_TEST_OPT_OBJS) -o $@ ${LIBS} ${LIBS_OPT}

${PERF_TEST_EXE}: $(PERF_TEST_OBJS) ${LIBS_OPT}
	$(CC) ${CFLAGS} ${LDFLAGS} $(PERF_TEST_OBJS) -o $@ ${LIBS} ${LIBS_OPT}

${LEGACY_PERF_TEST_EXE}: $(LEGACY_PERF_TEST_OBJS) ${LIBS_REF}
	$(CC) ${CFLAGS} ${LDFLAGS} $(LEGACY_PERF_TEST_OBJS) -o $@ ${LIBS} ${LIBS_REF}

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
