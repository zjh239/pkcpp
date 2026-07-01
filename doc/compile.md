## c++ version of Polikit

This is the C++ version of PAT (Polyhedral Analysis Tools). I try to use the c++20 standard and modules feature, so that the structure of the code can be quite similar to the Fortran version. Since modules feature is not fully supported by the compiler, some notes are put here.
#### Compile with CMake+Ninja

Since CMake 3.28+, it is possible to compile c++ code using the c++20 modules feature in a much easier way. To solve the problem of dependency chain of the code, a scanning is performed before the actual building. However, the Unix Makefiles generator does not support this scanning so far, and I use Ninja 1.12 in this case. The basic structure of `CMakeList.txt` in this case is

```cmake
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_SCAN_FOR_MODULES ON)
# first compile the modules
target_sources(my_app
    PRIVATE
	    FILE_SET modules TYPE CXX_MODULES FILES
		    ${MODULES}
)

# then compile the excutable file
target_sources(my_app PRIVATE main.cxx)
```

And the commands for building are

```bash
cmake -G Ninja ../src/
ninja
```

#### Compile std library

With modules feature, the standard libs can be included at the beginning of a module in following way:

```cpp
module;
#include <iostream>
export module my_mod;
```

If it is the main file, the last line is not needed. Instead of including the header of standard library, the std library need to be compiled to module before compiling our own code. As follows:

```bash
g++ -std=c++20 -fmodules-ts -x c++-system-header iostream
g++ -std=c++20 -fmodules-ts -x c++-system-header vector
```

Then they can be included as 
```cpp
import <iostream>;
import <vector>;
```
#### Compile code (with g++)

Suppose we have three files, `main.cpp` (the file that imports the module), `lib.cpp` (header or interface file of the module), and `lib.cxx` (implement or body of the module). With GNU compiler the compiling is

```bash
g++ -std=c++20 -fmodules-ts main.cpp lib.cpp lib.cxx
```

```bash
g++ -std=c++20 -fmodules-ts -c lib.cxx
g++ -std=c++20 -fmodules-ts -c main.cxx
g++ -std=c++20 -fmodules-ts *.o
```

Also, in some documents the `.cppm` suffix is used for the module files. Because it is not a recognizable suffix for GNU C++ compiler in default, the `-x c++` option must be used when compiling these files. If `.cxx` or `.cpp` suffix is used, there is no need to have `-x c++`.

```bash
g++ -std=c++20 -fmodules-ts main.cpp lib.cpp -x c++ lib.cppm
```

And, after creating `lib.o`, the conventional command also works.

```bash
g++ -std=c++20 -fmodules-ts main.cpp lib.o
```
#### Thoughts after the FORTRAN version

The old version is good, but I hope the data management and module dependency can be somehow improved.

I enjoyed the module feature of Fortran, that I could call functions without providing a long list of variables every time, but it turns out to be messy when there are more reusing of functions. Only when a certain module is called, the data structure is loaded.

In cpp, although I hope to use the module feature so that the code structure can be similar to fortran, they are still quite different. What I hope is to optimize at this moment is the I/O part and command line input part, so that all the computing modules can easily get the variables such as pair-wise cutoffs, and PBC conditions, and maybe the selection mask.

One more thing, in the cpp version I hope to implement parallel analysis from the beginning. The first part can be a neighbor finder.
