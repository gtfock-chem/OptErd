include ../../make.in

SRC =   erd__1111_csgto.f \
	erd__1111_def_blocks.f \
	erd__2d_atom_coefficients.f \
	erd__2d_atom_pq_integrals.f \
	erd__2d_coefficients.f \
	erd__2d_pq_integrals.f \
	erd__atom_int2d_to_e000.f \
	erd__atom_int2d_to_e0f0.f \
	erd__cartesian_norms.f \
	erd__csgto.f \
	erd__ctr_1st_half.f \
	erd__ctr_2nd_half_new.f \
	erd__ctr_2nd_half_update.f \
	erd__ctr_4index_block.f \
	erd__ctr_4index_reorder.f \
	erd__ctr_rs_expand.f \
	erd__ctr_tu_expand.f \
	erd__dsqmin_line_segments.f \
	erd__e0f0_def_blocks.f \
	erd__e0f0_pcgto_block.f \
	erd__gener_eri_batch.f \
	erd__hrr_matrix.f \
	erd__hrr_step.f \
	erd__hrr_transform.f \
	erd__int2d_to_e000.f \
	erd__int2d_to_e0f0.f \
	erd__map_ijkl_to_ikjl.f \
	erd__memory_1111_csgto.f \
	erd__memory_csgto.f \
	erd__memory_eri_batch.f \
	erd__move_ry.f \
	erd__normalize_cartesian.f \
	erd__pppp_pcgto_block.f \
	erd__rys_1_roots_weights.f \
	erd__rys_2_roots_weights.f \
	erd__rys_3_roots_weights.f \
	erd__rys_4_roots_weights.f \
	erd__rys_5_roots_weights.f \
	erd__rys_roots_weights.f \
	erd__rys_x_roots_weights.f \
	erd__set_abcd.f \
	erd__set_ij_kl_pairs.f \
	erd__spherical_transform.f \
	erd__sppp_pcgto_block.f \
	erd__sspp_pcgto_block.f \
	erd__sssp_pcgto_block.f \
	erd__ssss_pcgto_block.f \
	erd__transpose_batch.f \
	erd__xyz_to_ry_abcd.f \
	erd__xyz_to_ry_matrix.f \
	erd__prepare_ctr.F

SRC += erd_profile.c

ifeq "${arch}" "mic"
ERDLIB = liberd_mic.a
OBJS0 := $(addsuffix .o, $(SRC))
OBJS := $(addprefix mic/, $(OBJS0))
else
ERDLIB = liberd.a
OBJS := $(addsuffix .o, $(SRC))
endif

FFLAGS = -O3 -no-vec
FFLAGS += -reentrancy threaded -recursive
FFLAGS += ${OPTFLAGS}

CFLAGS = -O3 -std=gnu99 -openmp
CFLAGS += -Wall -w2 -Wunknown-pragmas -Wunused-variable
CFLAGS += ${OPTFLAGS}

all: ${ERDLIB}

${ERDLIB}: ${OBJS}
	${AR} $@ $^
	cp -r $@ ../lib/
	cp -r erd_profile.h ../include

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
