all: libs

libs:
	make -C external
	make -C libcint

scf:
	make -C pscf

test:
	make -C testprog

clean:
	make -C libcint clean
	make -C testprog clean

cleanall: clean
	make -C external clean
