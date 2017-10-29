# Nurgle

## How to

Put EcoCyc data files into `inputs` directory.

```
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=./local -DCMAKE_CXX_COMPILER=clang++ ..
$ make
$ make install
$ ./local/bin/main ./local/share/nurgle
```

## Installation

Packages, `clang-3.8` and `libc++-dev`, are required to build this.

When facing the error like:

```
/usr/include/c++/v1/cxxabi.h:21:10: fatal error: '__cxxabi_config.h' file not found
#include <__cxxabi_config.h>
         ^
         1 error generated.
```

install `libc++abi-dev` and do as follows: `sudo cp /usr/include/libcxxabi/__cxxabi_config.h /usr/include/c++/v1/.`.
