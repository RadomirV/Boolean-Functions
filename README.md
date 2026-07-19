# Boolean Functions

[![CMake CI](https://github.com/RadomirV/Boolean-Functions/actions/workflows/cmake.yml/badge.svg)](https://github.com/RadomirV/Boolean-Functions/actions/workflows/ci.yml)

A C++17 library for representing Boolean functions, computing cryptographic
properties, and experimenting with nonlinearity-improvement algorithms.

The project provides a `bf::BF` class, utility algorithms over bit vectors and
matrices over GF(2), a small example executable, and a GoogleTest test suite.

## Features

- Construction of Boolean functions from truth-table strings or generated data.
- Compact storage of truth tables in `uint32_t` blocks.
- Bit-level value access and mutation.
- Boolean-function metrics:
  - Hamming weight;
  - algebraic degree;
  - nonlinearity;
  - correlation immunity;
  - propagation and autocorrelation-related criteria.
- Mobius transform for algebraic normal form computations.
- Walsh-Hadamard transform.
- Autocorrelation transform.
- Basic linear algebra over GF(2):
  - Boolean matrix rank;
  - row-reduction-like transformations;
  - solving linear systems represented with bit masks.
- Experimental algorithms for finding value swaps that improve
  nonlinearity.

## Project Structure

```text
.
|-- CMakeLists.txt
|-- README.md
|-- src/
|   |-- BF.h
|   |-- BF.cpp
|   |-- BFTransforms.cpp
|   |-- BFMetrics.cpp
|   |-- BFGeneration.cpp
|   |-- BFImprove.cpp
|   |-- Utils.h
|   |-- Utils.cpp
|   `-- main.cpp
`-- test/
    `-- bf_tests.cpp
```

`BF.h` contains the public `bf::BF` interface. The implementation is split by
responsibility into core operations, transforms, metrics, generation, and
nonlinearity-improvement code. Utility algorithms live in `Utils.h/.cpp`.

## Requirements

- CMake 3.16 or newer.
- A C++17 compiler.
- Internet access during the first configure step, because GoogleTest is
  downloaded through CMake `FetchContent`.

On Windows, the project is currently configured and tested with Visual Studio
2022 generator settings.

## Build

Configure the project:

```powershell
cmake -S . -B build -DBUILD_TESTING=ON
```

Build all targets:

```powershell
cmake --build build --config Debug
```

Build only the example executable:

```powershell
cmake --build build --config Debug --target bf_main
```

Build only the tests:

```powershell
cmake --build build --config Debug --target bf_tests
```

## Run Tests

```powershell
ctest --test-dir build -C Debug --output-on-failure
```

The test suite covers bit utilities, GF(2) matrix routines, core `bf::BF`
operations, Mobius transform, Walsh-Hadamard transform, and nonlinearity-related
checks.

## Continuous Integration

GitHub Actions builds the project on Linux and Windows using CMake. The workflow
configures the project, builds `bf_tests` and `bf_main`, and runs the test suite
with `ctest`.


## Example

```cpp
#include "BF.h"

#include <iostream>
#include <vector>

int main()
{
    bf::BF f("0001");

    std::cout << "Variables: " << f.getNumberOfVariables() << '\n';
    std::cout << "Weight: " << f.weight() << '\n';

    std::vector<int> wht = f.WHTransform();
    std::cout << "Nonlinearity: " << f.nonlinearity(wht) << '\n';

    bf::BF anf = f.mobiusTransform();
    std::cout << "Degree: " << anf.degree() << '\n';
}
```
