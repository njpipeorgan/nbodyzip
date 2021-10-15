# nbodyzip

A lossy compressor of dark matter-only N-body simulation data

## Installation and usage

To compile the source code,

```bash
g++ -std=c++17 -o nbodyzip main.cpp
```

To use nbodyzip,

```bash
./nbodyzip [input_format] [output_format] [input_file] [output_file]
```

The input/output format can be one of `tipsy`, `posvel`, and `zip`.

For example, to compress a pkdgrav3 snapshot, run

```bash
./nbodyzip tipsy zip snapshot.xxxxx snapshot.xxxxx.nz
```

To decompress, run

```bash
./nbodyzip zip tipsy snapshot.xxxxx.nz snapshot.xxxxx.2
```

## Input/output formats

### `tipsy`

nbodyzip supports `tipsy` files that contains only dark matter particles and the number of particles along each dimension should be a power of two.

### `posvel`

is a binary format (system-dependent endianness) defined in nbodyzip. It contains a 32-byte header, an array of particle positions, and an array of particle velocities.

- Header: scale factor (fp64), number of particles (ui64), particle mass (fp64), and softening length (fp64, in units of mean particle distance);
- Positions: *x*<sub>1</sub>, *y*<sub>1</sub>, *z*<sub>1</sub>, ..., *x*<sub>*n*</sub>, *y*<sub>*n*</sub>, *z*<sub>*n*</sub> (ui32 numbers representing coordinates between 0 and 1)
- Velocites: *v*<sub>*x*1</sub>, *v*<sub>*y*1</sub>, *v*<sub>*z*1</sub>, ..., *v*<sub>*xn*</sub>, *v*<sub>*yn*</sub>, *v*<sub>*zn*</sub> (fp32)

### `zip`

is a binary format (system-dependent endianness) defined in nbodyzip that efficiently stores simulation data. For a 1024<sup>3</sup> simulation, the absolute errors of *x*, *y*, and *z* coordinates are less than 2<sup>-22</sup>, 2<sup>-22</sup>, and  2<sup>-21</sup>; the absolute error of *v*<sub>*i*</sub> is 1.5&times;10<sup>-5</sup>(1+(*v*<sub>*i*</sub>/*v*<sub>scale</sub>)), with *v*<sub>scale</sub> being the rms of *v*<sub>*i*</sub>.
