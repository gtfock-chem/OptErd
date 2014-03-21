all: libs

libs:
	make -C external
	make -C libcint

clean:
	make -C libcint clean

cleanall: clean
	make -C external clean
