include ../make.in

EXEC = test_fock

CFLAGS = -O3 -Wall -openmp -w2 -std=gnu99
CFLAGS += -Wunknown-pragmas -Wunused-variable -Wunused-function
CFLAGS += ${OPTFLAGS} #-vec-report2

INC = -I./  \
      -I../include/ \
      -I../external/Yeppp/include
LDFLAGS = -L./ -L../lib/

SRC = $(wildcard *.c)
OBJS := $(addsuffix .o, $(basename $(SRC)))
ifeq "${arch}" "mic"
LIBS = -lcint_mic -lerd_mic -loed_mic
else
LIBS = ../lib/libcint.a ../lib/liberd.a ../lib/liboed.a
endif

all: ${EXEC}

${EXEC}: $(OBJS)
	$(CC) ${CFLAGS} ${LDFLAGS} $(OBJS) -o $@ ${LIBS}

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

clean:
	rm -f *.o *.s *.d *~ ${EXEC}
