# Boolean Functions

A small C++17 project for representing Boolean functions major cryptographic properties and methods for improving nonlinearity.

The project contains a `BF` class, utility algorithms over bit vectors and
matrices over GF(2), and a GoogleTest-based tests.

## Features

- Construction of Boolean functions from truth-table strings.
- Bit-level access and mutation of function values.
- Boolean function metrics:
  - Hamming weight;
  - algebraic degree;
  - nonlinearity;
  - correlation-related criteria;
  - propagation/autocorrelation-related criteria.
- Mobius transform for algebraic normal form computations.
- Walsh-Hadamard transform.
- Autocorrelation transform.
- Basic linear algebra over GF(2):
  - Boolean matrix rank;
  - row-reduction-like transformations;
  - solving linear systems represented with bit masks.
- Experimental algorithms for finding value swaps that may improve
  nonlinearity.

Some names and module boundaries are still being improved. See
[Roadmap](#roadmap).

## Project Structure

```text
.
|-- CMakeLists.txt
|-- README.md
|-- src/
|   |-- BF.h
|   |-- BF.cpp
|   `-- main.cpp
`-- test/
    `-- bf_tests.cpp
```

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

The test suite covers bit utilities, GF(2) matrix routines, core `BF`
operations, Mobius transform, Walsh-Hadamard transform, and nonlinearity-related
checks.

## Example

```cpp
#include "BF.h"

#include <iostream>
#include <vector>

int main()
{
    BF f("0001");

    std::cout << "Variables: " << f.getNumberOfVariables() << '\n';
    std::cout << "Weight: " << BF::weight(f) << '\n';

    std::vector<int> wht = BF::WHTransform(f);
    std::cout << "nonlinearity: " << BF::nonlinearity(f, wht) << '\n';
}
```

## Notes

Boolean function values are stored compactly in `uint32_t` blocks. Many
algorithms operate directly on bit masks, which keeps the implementation close
to the mathematical representation used in Boolean-function research.

Several APIs are intentionally still simple and experimental. The current goal
is to keep the implementation testable while gradually improving naming,
module boundaries, and error handling.

## Roadmap

- Split utility algorithms out of `BF.h` into dedicated modules:
  - bit utilities;
  - GF(2) linear algebra;
  - Boolean-function transforms;
  - Boolean-function metrics.
- Rename remaining abbreviated APIs to clearer names, for example
  `WHTransform`, `autoCor`, `correlationImmunity`, `propagationCriteria`, and `cnF`.
- Replace console diagnostics in library code with exceptions or explicit error
  handling.
- Add `.clang-format`.
- Add GitHub Actions for CMake build and tests.
- Add more examples and documentation for the nonlinearity-improvement
  algorithms.

## License

See [LICENSE.txt](LICENSE.txt).
