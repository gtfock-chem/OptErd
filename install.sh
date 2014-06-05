python configure.py $1
mkdir -p ${prefix}
rm -rf ${prefix}/lib ${prefix}/include
mkdir -p ${prefix}/lib
cp -f lib/$1/liberd-opt.a ${prefix}/lib/liberd.a
cp -f lib/$1/liboed-opt.a ${prefix}/lib/liboed.a
cp -f lib/$1/libcint-opt.a ${prefix}/lib/libcint.a
cp -rf include/ ${prefix}
