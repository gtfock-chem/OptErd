include ../../make.in

SRC = $(wildcard *.c)

ifeq "${arch}" "mic"
ERDLIB = liberd_mic.a
OBJS0 := $(addsuffix .o, $(SRC))
OBJS := $(addprefix mic/, $(OBJS0))
else
ERDLIB = liberd.a
OBJS := $(addsuffix .o, $(SRC))
endif

FFLAGS = -O3
FFLAGS += ${OPTFLAGS}

CFLAGS = -O3 -std=gnu99 -I../Yeppp/include -fopenmp
CFLAGS += -Wall -Wunused-variable
CFLAGS += ${OPTFLAGS}
CFLAGS += -D__ALIGNLEN__=${alignlen}

all: ${ERDLIB}
	cp erd_profile.h ${prefix}/include

${ERDLIB}: ${OBJS}
	${AR} $@ $^
	cp -r *.a ${prefix}/lib/

%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

%.f.o: %.f
	${FC} ${FFLAGS} -c $< -o $@

%.F.o: %.F
	${FC} ${FFLAGS} -c $< -o $@

mic/%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

mic/%.f.o: %.f
	${FC} ${FFLAGS} -c $< -o $@

mic/%.F.o: %.F
	${FC} ${FFLAGS} -c $< -o $@

clean:
	rm -f *.o mic/*.o *~ *.a
