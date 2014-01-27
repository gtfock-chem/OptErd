#!/usr/bin/python

from __future__ import print_function

import optparse
import os
import sys
import glob

root_dir = os.path.dirname(__file__)

erd_f90_sources = [os.path.relpath(path, root_dir) for path in glob.glob(os.path.join(root_dir, 'external', 'erd', '*.f'))]
erd_f77_sources = [os.path.relpath(path, root_dir) for path in glob.glob(os.path.join(root_dir, 'external', 'erd', '*.F'))]
erd_ref_sources = erd_f90_sources + erd_f77_sources + [os.path.join("external", "erd", "erd_profile.c")]
erd_opt_sources = [os.path.join("external", "erd", filename) for filename in [
	"erd__1111_csgto.c", "erd__1111_def_blocks.c", "erd__2d_coefficients.c", "erd__2d_pq_integrals.c", "erd__cartesian_norms.c", "erd__csgto.c",
	"erd__ctr_4index_block.c", "erd__dsqmin_line_segments.c", "erd__e0f0_def_blocks.c", "erd__e0f0_pcgto_block.c", "erd__gener_eri_batch.f",
	"erd__hrr_matrix.c", "erd__hrr_step.c", "erd__hrr_transform.c", "erd__int2d_to_e000.c", "erd__int2d_to_e0f0.c",
	"erd__memory_1111_csgto.c", "erd__memory_csgto.c", "erd__memory_eri_batch.c", "erd__move_ry.c", "erd__normalize_cartesian.c",
	"erd__pppp_pcgto_block.c", "erd__rys_1_roots_weights.f", "erd__rys_2_roots_weights.f", "erd__rys_3_roots_weights.f", "erd__rys_4_roots_weights.f",
	"erd__rys_5_roots_weights.f", "erd__rys_roots_weights.f", "erd__rys_x_roots_weights.f", "erd__set_abcd.c", "erd__set_ij_kl_pairs.c",
	"erd__spherical_transform.c", "erd__sppp_pcgto_block.c", "erd__sspp_pcgto_block.c", "erd__sssp_pcgto_block.c", "erd__ssss_pcgto_block.c",
	"erd__xyz_to_ry_abcd.c", "erd__xyz_to_ry_matrix.c", "erd__prepare_ctr.c", "erd__boys_table.c", "erd_profile.c"
]] + erd_f77_sources
erd_c_sources = filter(lambda source_file: source_file.endswith('.c'), erd_opt_sources)

oed_f90_sources = [os.path.relpath(path, root_dir) for path in glob.glob(os.path.join(root_dir, 'external', 'oed', '*.f'))]
oed_f77_sources = [os.path.relpath(path, root_dir) for path in glob.glob(os.path.join(root_dir, 'external', 'oed', '*.F'))]

cint_sources = [os.path.relpath(path, root_dir) for path in glob.glob(os.path.join(root_dir, 'libcint', '*.c'))]

tab = '  '

with open('build.ninja', 'w') as makefile:
	print('FC_NHM = ifort -m64 -xSSE4.2', file = makefile)
	print('FC_SNB = ifort -m64 -xAVX', file = makefile)
	print('FC_MIC = ifort -mmic', file = makefile)
	print('FFLAGS = -O3 -g -reentrancy threaded -recursive', file = makefile)
	print('CC_NHM = icc -m64 -xSSE4.2', file = makefile)
	print('CC_SNB = icc -m64 -xAVX', file = makefile)
	print('CC_MIC = icc -mmic -no-opt-prefetch', file = makefile)
	print('CFLAGS = -O3 -g -std=gnu99 -D__ERD_PROFILE__ -Iexternal/Yeppp/include -w -Wunknown-pragmas -Wno-unused-variable', file = makefile)
	print('LDFLAGS = -static-intel -lifcore -openmp', file = makefile)
	print('AR = xiar', file = makefile)

	print('rule COMPILE_F90', file = makefile)
	print(tab + 'command = $FC $FFLAGS -o $out -c $in', file = makefile)
	print(tab + 'description = F90[$ARCH] $in', file = makefile)

	print('rule COMPILE_F77', file = makefile)
	print(tab + 'command = $FC $FFLAGS -o $out -c $in', file = makefile)
	print(tab + 'description = F77[$ARCH] $in', file = makefile)

	print('rule COMPILE_C', file = makefile)
	print(tab + 'depfile = $DEP_FILE', file = makefile)
	print(tab + 'command = $CC $CFLAGS -MMD -MT $out -MF $DEP_FILE -o $out -c $in', file = makefile)
	print(tab + 'description = CC[$ARCH] $in', file = makefile)

	print('rule LINK', file = makefile)
	print(tab + 'command = $CC $CFLAGS -o $out $in $LDFLAGS', file = makefile)
	print(tab + 'description = CCLD[$ARCH] $in', file = makefile)

	print('rule CREATE_STATIC_LIBRARY', file = makefile)
	print(tab + 'command = $AR rcs $out $in', file = makefile)
	print(tab + 'description = AR[$ARCH] $out', file = makefile)

	for arch in ['nhm', 'snb', 'mic']:
		for source_file in erd_f90_sources:
			object_file = os.path.join(os.path.dirname(source_file), arch, source_file + '.o')
			print('build %s : COMPILE_F90 %s' % (object_file, source_file), file = makefile)
			print(tab + 'FC = $FC_%s' % arch.upper(), file = makefile)
			print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		for source_file in erd_f77_sources:
			object_file = os.path.join(os.path.dirname(source_file), arch, source_file + '.o')
			print('build %s : COMPILE_F77 %s' % (object_file, source_file), file = makefile)
			print(tab + 'FC = $FC_%s' % arch.upper(), file = makefile)
			print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		for source_file in erd_c_sources:
			object_file = os.path.join(os.path.dirname(source_file), arch, source_file + '.o')
			dep_file = os.path.join(os.path.dirname(source_file), arch, source_file + '.d')
			print('build %s : COMPILE_C %s' % (object_file, source_file), file = makefile)
			print(tab + 'DEP_FILE = ' + dep_file, file = makefile)
			print(tab + 'CC = $CC_%s' % arch.upper(), file = makefile)
			print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		erd_ref_objects = [os.path.join(os.path.dirname(source_file), arch, source_file + '.o') for source_file in erd_ref_sources]
		erd_opt_objects = [os.path.join(os.path.dirname(source_file), arch, source_file + '.o') for source_file in erd_opt_sources]

		print('build lib/' + arch + '/liberd_ref.a : CREATE_STATIC_LIBRARY ' + ' '.join(erd_ref_objects), file = makefile)
		print(tab + 'ARCH = %s' % arch.upper(), file = makefile)
		print('build lib/' + arch + '/liberd_opt.a : CREATE_STATIC_LIBRARY ' + ' '.join(erd_opt_objects), file = makefile)
		print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		for source_file in oed_f90_sources:
			object_file = os.path.join(os.path.dirname(source_file), arch, source_file + '.o')
			print('build %s : COMPILE_F90 %s' % (object_file, source_file), file = makefile)
			print(tab + 'FC = $FC_%s' % arch.upper(), file = makefile)
			print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		for source_file in oed_f77_sources:
			object_file = os.path.join(os.path.dirname(source_file), arch, source_file + '.o')
			print('build %s : COMPILE_F77 %s' % (object_file, source_file), file = makefile)
			print(tab + 'FC = $FC_%s' % arch.upper(), file = makefile)
			print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		oed_sources = oed_f90_sources + oed_f77_sources
		oed_objects = [os.path.join(os.path.dirname(source_file), arch, source_file + '.o') for source_file in oed_sources]
		print('build lib/' + arch + '/liboed.a : CREATE_STATIC_LIBRARY ' + ' '.join(oed_objects), file = makefile)
		print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		for source_file in cint_sources:
			object_file = os.path.join(os.path.dirname(source_file), arch, source_file + '.o')
			dep_file = os.path.join(os.path.dirname(source_file), arch, source_file + '.d')
			print('build %s : COMPILE_C %s' % (object_file, source_file), file = makefile)
			print(tab + 'DEP_FILE = ' + dep_file, file = makefile)
			print(tab + 'CC = $CC_%s' % arch.upper(), file = makefile)
			print(tab + 'ARCH = %s' % arch.upper(), file = makefile)
		cint_objects = [os.path.join(os.path.dirname(source_file), arch, source_file + '.o') for source_file in cint_sources]
		print('build lib/' + arch + '/libcint.a : CREATE_STATIC_LIBRARY ' + ' '.join(cint_objects), file = makefile)
		print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		screening_object_file = 'testprog/%s/screening.c.o' % arch
		print('build %s : COMPILE_C testprog/screening.c' % screening_object_file, file = makefile)
		print(tab + 'DEP_FILE = testprog/%s/screening.c.d' % arch, file = makefile)
		print(tab + 'CC = $CC_%s' % arch.upper(), file = makefile)
		print(tab + 'CFLAGS = $CFLAGS -Iinclude -openmp', file = makefile)
		print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		for variant in ['opt', 'ref']:
			object_file = 'testprog/%s/testCInt.c.%s.o' % (arch, variant)
			print('build %s : COMPILE_C testprog/testCInt.c' % object_file, file = makefile)
			print(tab + 'DEP_FILE = %s.d' % object_file, file = makefile)
			print(tab + 'CC = $CC_%s' % arch.upper(), file = makefile)
			print(tab + 'CFLAGS = $CFLAGS -Iexternal/erd -Iinclude -D' + {'opt': 'OPTERD_TEST_OPTIMIZED', 'ref': 'OPTERD_TEST_REFERENCE'}[variant], file = makefile)
			print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

			binary_file = 'testprog/%s/Test.%s' % (arch, variant.title())
			libs = " ".join(['lib/' + arch + '/lib' + lib + '.a' for lib in ['cint', 'oed', 'erd_' + variant]])
			print('build %s : LINK %s %s %s' % (binary_file, object_file, screening_object_file, libs), file = makefile)
			print(tab + 'CC = $CC_%s' % arch.upper(), file = makefile)
			print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		object_file = 'testprog/%s/testCInt2.c.o' % arch
		print('build %s : COMPILE_C testprog/testCInt2.c' % object_file, file = makefile)
		print(tab + 'DEP_FILE = testprog/%s/testCInt2.c.d' % arch, file = makefile)
		print(tab + 'CC = $CC_%s' % arch.upper(), file = makefile)
		print(tab + 'CFLAGS = $CFLAGS -Iexternal/erd -Iinclude -openmp', file = makefile)
		print(tab + 'ARCH = %s' % arch.upper(), file = makefile)

		binary_file = 'testprog/%s/Test.Perf' % arch
		libs = " ".join(['lib/' + arch + '/lib' + lib + '.a' for lib in ['cint', 'oed', 'erd_opt']])
		print('build %s : LINK %s %s %s' % (binary_file, object_file, screening_object_file, libs), file = makefile)
		print(tab + 'CC = $CC_%s' % arch.upper(), file = makefile)
		print(tab + 'ARCH = %s' % arch.upper(), file = makefile)
