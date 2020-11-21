swig -python -py3 raxml.i
gcc -c -fpic raxml_wrap.c -I /usr/include/python3.5 -I src -o _raxml.o
gcc -shared _raxml.o -o _raxml.so
