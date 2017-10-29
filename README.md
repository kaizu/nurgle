# Nurgle

> "Buboes, phlegm, blood and guts! Boils, bogeys, rot and pus! Blisters, fevers, weeping sores! From your wounds the fester pours."
> 
> -- <cite>The Chant of Nurgle</cite>

## How to

Packages, `clang-3.8` and `libc++-dev`, and Boost headers are required to build this.

Put EcoCyc data files into `inputs/21.1` directory, and build as follows.

```
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=./local -DCMAKE_CXX_COMPILER=clang++ ..
$ make
$ make install
$ ./local/bin/main ./local/share/nurgle
```

## Troubleshooting

When facing the error like:

```
/usr/include/c++/v1/cxxabi.h:21:10: fatal error: '__cxxabi_config.h' file not found
#include <__cxxabi_config.h>
         ^
         1 error generated.
```

install `libc++abi-dev` and try: `sudo cp /usr/include/libcxxabi/__cxxabi_config.h /usr/include/c++/v1/.`.
