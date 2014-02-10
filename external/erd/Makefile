include ../../make.in

SRC   = erd__1111_csgto.c \
	erd__1111_def_blocks.c \
	erd__2d_coefficients.c \
	erd__2d_pq_integrals.c \
	erd__boys_table.c \
	erd__jacobi_table.c \
	erd__cartesian_norms.c \
	erd__csgto.c \
	erd__dsqmin_line_segments.c \
	erd__e0f0_def_blocks.c \
	erd__e0f0_pcgto_block.c \
	erd__hrr_matrix.c \
	erd__hrr_step.c \
	erd__hrr_transform.c \
	erd__int2d_to_e000.c \
	erd__int2d_to_e0f0.c \
	erd__memory_1111_csgto.c \
	erd__memory_csgto.c \
	erd__move_ry.c \
	erd__normalize_cartesian.c \
	erd__pppp_pcgto_block.c \
	erd__rys_1_roots_weights.c \
	erd__rys_2_roots_weights.c \
	erd__rys_3_roots_weights.c \
	erd__rys_4_roots_weights.c \
	erd__rys_5_roots_weights.c \
	erd__rys_roots_weights.c \
	erd__rys_x_roots_weights.c \
	erd__set_abcd.c \
	erd__set_ij_kl_pairs.c \
	erd__spherical_transform.c \
	erd__sppp_pcgto_block.c \
	erd__sspp_pcgto_block.c \
	erd__sssp_pcgto_block.c \
	erd__ssss_pcgto_block.c \
	erd__xyz_to_ry_abcd.c \
	erd__xyz_to_ry_matrix.c 

SRC += erd_profile.c

ifeq "${arch}" "mic"
ERDLIB = liberd_mic.a
OBJS0 := $(addsuffix .o, $(SRC))
OBJS := $(addprefix mic/, $(OBJS0))
else
ERDLIB = liberd.a
OBJS := $(addsuffix .o, $(SRC))
endif

FFLAGS = -O3
FFLAGS += -reentrancy threaded -recursive
FFLAGS += ${OPTFLAGS}

CFLAGS = -O3 -std=gnu99 -I../Yeppp/include -openmp
CFLAGS += -Wall -w2 -Wunknown-pragmas -Wunused-variable
CFLAGS += ${OPTFLAGS} -restrict
CFLAGS += -D__ALIGNLEN__=${alignlen}
CFLAGS += -vec-report=2

all: ${ERDLIB}
	cp erd_profile.h ../../include

${ERDLIB}: ${OBJS}
	${AR} $@ $^
	cp -r *.a ../../lib/

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
