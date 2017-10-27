# Nurgle

## How to

Put EcoCyc data files into `inputs` directory.

```
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=./local -DCMAKE_CXX_COMPILER=clang++ ..
$ make
$ make install
$ ./local/bin/main
```
