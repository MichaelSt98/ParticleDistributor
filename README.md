# ParticleDistributor

Create **Particle Distribution(s)** for N-Body (and prospectively SPH) simulation.

## Dependencies

The following dependencies need to be installed to build this program.

### HDF5

| library         | licence           | usage             | link               |
| --------------- | ----------------- | ----------------- | ------------------ |
| HDF5            | HDF5 License (BSD-Style) | HDF5 for I/O operations | [hdf5group.org](https://www.hdfgroup.org/solutions/hdf5/) |

Can be installed via a package manager e.g.

```
$ sudo apt-get install libhdf5-serial-dev
```

**or** build from source:

* `<hdf5 version>` e.g. `1.12.2`

```
$ wget https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_2/source/hdf5-<hdf5 version>.tar.gz
$ tar zxvf hdf5-<hdf5 version>.tar.gz
$ cd hdf5-<hdf5 version>
$ ./configure --prefix=<install-directory>
$ make			# build the library
$ make check	# verify the correctness
$ make install
```
* taken from [realease_docs](https://github.com/HDFGroup/hdf5/tree/develop/release_docs)

The `<install-directory>` has to be set in the Makefile via the variable `HDF5DIR`.

### Boost
Can be installed via a package manager e.g.

```
$ sudo apt-get install libboost-all-dev
```

**or** build from source:


| library         | licence           | usage             | link               |
| --------------- | ----------------- | ----------------- | ------------------ |
| Boost           | Boost Software License 1.0 | config file parsing | [boost.org](https://www.boost.org/) |

* `<boost version>` e.g. `1_78_0`

```
$ wget https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_<boost version>.tar.gz
$ tar zxvf boost_<boost version>.tar.gz
$ cd boost_<boost version>
$ ./bootstrap.sh --with-libraries=all
$ ./b2
$ ./b2 install --prefix=<install-directory>
```

### HighFive

| library         | licence           | usage             | link               |
| --------------- | ----------------- | ----------------- | ------------------ |
| HighFive        | Boost Software License 1.0 | C++ wrapper for HDF5 | [github.com/BlueBrain/HighFive](https://github.com/BlueBrain/HighFive) |

* `git clone https://github.com/BlueBrain/HighFive.git`

Copy the folder `include/highfive` to `HEADERONLYDIR` defined in the Makefile

### Cxxopts

| library         | licence           | usage             | link               |
| --------------- | ----------------- | ----------------- | ------------------ |
| cxxopts         | MIT license | command line argument parsing | [github.com/jarro2783/cxxopts](https://github.com/jarro2783/cxxopts) |

* `git clone https://github.com/jarro2783/cxxopts.git`

Copy the file `include/cxxopts.hpp` to `HEADERONLYDIR` defined in the Makefile

## Build

Update the following variables in the Makefile as described in the dependencies section:

* `HEADERONLYDIR`
* `HDF5DIR`

Then you can build the program with `make`.

## Getting started

**Note:** If you installed HDF5 locally you may need to update the environment variable `LD_LIBRARY_PATH` to include your HDF5 installation. This can be done in a bash shell with:

```
export LD_LIBRARY_PATH=$(HDF5DIR)/lib:$(LD_LIBRARY_PATH)
```
The program's help can then be showed via `bin/runner -h` which gives:

```
Generating HDF5 file with initial particle distribution for N-Body and SPH simulation.
Usage:
  bin/runner [OPTION...]

  -N, --N-particles arg       Number of particles (default: 1000000)
  -s, --seed arg              Use given random seed (default: 0)
  -f, --filename arg          File name
  -d, --distributionType arg  ([0]: plummer, [1]: single_spiral, [2]: kepler, [3]: 
                              multiple_plummer, [4]: kepler_torus) (default: 0)
  -h, --help                  Show this help
```

## Distributors

**Currently supported distributors:**

* [Plummer](src/distributors/Plummer.cpp)
* [SingleSpiral](src/distributors/SingleSpiral.cpp)
* [Kepler](src/distributors/Kepler.cpp)
* [MultiplePlummer](src/MultiplePlummer.cpp)
* [KeplerTorus](src/KeplerTorus.cpp)

To change the configuration of the distributor edit the corrensponding configuration file found in the [config folder](config/).

### Plummer

**Plummer model**

A plummer sphere with 10000 particles and a random seed of 12345 can be generated with

```
$ mkdir output
$ bin/runner -N 10000 -d 0 -f pl -s 12345
```
which produces an output file at `output/plN10000seed12345.h5`.

### SingleSpiral

**Two rigid rotating spheres evolving into some kind of spiral**


### Kepler

**Central star with small bodies orbiting with kepler velocity**


### MultiplePlummer

**Multiple Plummer spheres with their centers randomly distributed in a box** 

### KeplerTorus

**Central star with small bodies distributed on a torus orbiting with kepler velocity**