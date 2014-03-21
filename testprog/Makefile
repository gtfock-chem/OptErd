include ../make.in

CFLAGS = -O3 -Wall -openmp -w2 -std=gnu99
CFLAGS += -Wunknown-pragmas -Wunused-variable -Wunused-function
CFLAGS += ${OPTFLAGS}

INC = -I./  \
      -I../include/ \
      -I../external/Yeppp/include

LDFLAGS = -L./ -lifcore

ifeq "${arch}" "mic"
UNIT_TEST_EXE = Test.Opt.MIC
PERF_TEST_EXE = Test.Perf.MIC
UNIT_TEST_OBJS := mic/testCInt.o mic/screening.o
PERF_TEST_OBJS := mic/testCInt2.o mic/screening.o
LIBS = ../lib/libcint_mic.a ../lib/liboed_mic.a ../lib/liberd_mic.a
else
UNIT_TEST_EXE = Test.Opt
PERF_TEST_EXE = Test.Perf
UNIT_TEST_OBJS := testCInt.o screening.o
PERF_TEST_OBJS := testCInt2.o screening.o
LEGACY_PERF_TEST_OBJS := testCInt2.o screening.o
LIBS = ../lib/libcint.a ../lib/liboed.a ../lib/liberd.a
endif

all: ${UNIT_TEST_EXE} ${PERF_TEST_EXE}

${UNIT_TEST_EXE}: $(UNIT_TEST_OBJS) ${LIBS}
	$(CC) ${CFLAGS} ${LDFLAGS} $(UNIT_TEST_OBJS) -o $@ ${LIBS}
	cp -r $@ ../bin/

${PERF_TEST_EXE}: $(PERF_TEST_OBJS) ${LIBS}
	$(CC) ${CFLAGS} ${LDFLAGS} $(PERF_TEST_OBJS) -o $@ ${LIBS}
	cp -r $@ ../bin/

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

mic/%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

clean:
	rm -f *.o *.s *.d *~ mic/*.o testCInt*.o Test.*
