# Install MOOSE with Pytorch C++ API (libtorch)

The way one enables pytorch capabilities in MOOSE depends on
the operating system (Linux or Mac) and if we use HPC or just a local workstation.

!alert! note
Before we review the main approaches, it is important to emphasize that
linking MOOSE with libtorch on Linux machines is not support if the compiler stack has been built
using a `libc` version below 2.27 (for `libtorch v 1.8+`)
or 2.23 (for `libtorch v1.4-1.8`). Furthermore, we do not support `libtorch` versions below
v1.4. To check your currently used libc version use the following command:

```bash
ldd --version
```

!alert-end!

## Installation on a local Mac workstation

The first step of this process is to set up the a suitable conda environment. For this follow
the instructions [here](https://mooseframework.inl.gov/getting_started/installation/conda.html).
Next, activate the conda environment (assume `moose-torch`), clone the repository
 and navigate to the base folder of the
stochastic tools module (let's say it is in `~/projects`):

```bash
conda activate moose-torch
git clone https://github.com/grmnptr/moose.git
cd moose
git checkout torch-compile-19571
cd ~/projects/moose/modules/stochastic_tools
```

As a next step, we download a precompiled version of `libtorch` which is suitable
for our system. This is carried out by a bash script as:

```bash
./scripts/setup_libtorch.sh
```

!alert! note
The desired version of libtorch can be set by the following
argument:

```bash
./scripts/setup_libtorch.sh --version=1.8
```

!alert-end!

which creates a `libtorch` folder in the stochastic tools root directory. This folder
will contain the necessary header files and shared object files which can be used
for dynamic linking.

The last step is to compile our MOOSE with libtorch. For this, a non-unity build
is used (which is considerably slower than the default unity build) in the
following manner:

```bash
MOOSE_UNITY=false make -j 8
```

## Installation on a local Linux workstation

The difference compared to Mac machines is that the moose conda packages (compiler stack) have been
built using `libc` v2.12 meaning that it cannot be used to set up the framework with
libtorch. For this reason, we need to install PETSc and Libmesh manually, assuming that
the default compiler stack on the system uses libc 2.27+ (most modern distributions,
such as Ubuntu actually do have suitable compilers). For this first, clone the
repository, then use the corresponding scripts as:

```bash
git clone https://github.com/grmnptr/moose.git
cd moose
git checkout torch-compile-19571
./scripts/update_and_rebuild_petsc.sh
./scripts/update_and_rebuild_libmesh.sh
```

From this point on, the procedure is more or less the same as the installation
on a Mac machine. We use the same script to download a suitable
precompiled version of `libtorch`:

```bash
cd modules/stochastic_tools
./scripts/setup_libtorch.sh
```

which creates a `libtorch` folder in the stochastic tools root directory. This folder
will contain the necessary header files and shared object files which can be used
for dynamic linking.

Now, we can move on compile MOOSE with libtorch. Again, a non-unity build
is used to make sure that namespace conflicts between `libmesh` and `libtorch` are
avoided:

```bash
MOOSE_UNITY=false make -j 8
```

## Installation on a HPC machine

The installation process is very similar to the one on a local Linux machine.
However, in some cases, we may have to build our own compiler stack/mpicc.
For this, follow the instructions
[here](https://mooseframework.inl.gov/getting_started/installation/manual_installation_gcc.html).
Once this is done and we ensured that out compilers were built using glibc 2.7+,
we can proceed to clone moose, built `petsc` and `libmesh`, download `libtorch` and
built moose with `libtorch`:

```bash
git clone https://github.com/grmnptr/moose.git
cd moose
git checkout torch-compile-19571
./scripts/update_and_rebuild_petsc.sh
./scripts/update_and_rebuild_libmesh.sh
cd modules/stochastic_tools
./scripts/setup_libtorch.sh
MOOSE_UNITY=false make -j 8
```

## Testing the installation

A basic neural network trainer class has been implemented in the stochastic tools
module for testing purposes. If MOOSE compiled without error messages, we can
navigate to the testing folder within the stochastic tools root
and use the following command to train a simple neural network:

```bash
cd test/tests/surrogates/basic_nn/
./../../../../stochastic_tools-opt -i train-nn.i --allow-test-objects
```

The expected output should look like this:

```bash
Framework Information:
MOOSE Version:           git commit d4f7f666e3 on 2021-12-16
LibMesh Version:         
PETSc Version:           3.15.1
SLEPc Version:           3.15.1
Current Time:            Tue Dec 21 11:09:04 2021
Executable Timestamp:    Tue Dec 21 11:07:46 2021

Parallelism:
  Num Processors:          1
  Num Threads:             1

Mesh:
  Parallel Type:           replicated
  Mesh Dimension:          1
  Spatial Dimension:       1
  Nodes:                   2
  Elems:                   1
  Num Subdomains:          1

Execution Information:
  Executioner:             Steady
  Solver Mode:             Preconditioned JFNK

 Solve Skipped!
Epoch: 10 | Loss: 59.7538
Epoch: 20 | Loss: 56.9033
Epoch: 30 | Loss: 55.8542
Epoch: 40 | Loss: 54.5681
Epoch: 50 | Loss: 52.8059
Epoch: 60 | Loss: 50.5726
Epoch: 70 | Loss: 47.6894
Epoch: 80 | Loss: 44.2285
Epoch: 90 | Loss: 40.4647
Epoch: 100 | Loss: 36.7318
Epoch: 110 | Loss: 33.1301
...
...
```