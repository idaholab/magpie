# Installation

## Install and build MOOSE

Refer to the [MOOSE Getting Started](http://mooseframework.org/getting-started/)
pages to set up MOOSE System Environment and check out MOOSE from its git
repository. Magpie expects moose to either be checked out side by side

```
./projects/moose
./projects/magpie
```

or the `MOOSE_DIR` environment variable to point to the MOOSE directory.

## Install prerequisite packages

Magpie needs GSL (GNU Scientific Library) and the fast Fourier transform library
FFTW3 installed for a feature complete build. The configuration of these
external packages is performed using
[pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/).
Installation of these dependencies is system specific:

### Conda (cross platform)

This is now the **default** way to install MOOSE/Magpie dependencies.

```
conda install pkg-config gsl fftw
```

### macOS (Homebrew)

Prerequisite package installation on Mac is easiest using the [Homebrew](https://brew.sh/)
package manager:

```
brew install pkg-config fftw gsl
```

### Linux (Ubuntu / Debian)

Use the native package manager of your Linux distribution. On Debian based
systems, such as Ubuntu, use *apt* as follows:

```
sudo apt install pkg-config libgsl-dev libfftw3-dev
```

## Install Magpie

Magpie is hosted on [GitHub](https://github.com/idaholab/magpie) and can be
cloned directly from there using [Git](https://git-scm.com/). We recommend
creating a directory named `projects` to put all of your MOOSE related work
(such as Magpie) in which leads to the following commands (from your home
directory):

```bash
mkdir ~/projects
cd ~/projects
git clone https://github.com/idaholab/magpie.git
cd ~/projects/magpie
git checkout master
```

!alert note
The "master" branch of Magpie is the "stable" branch that will only be updated
after all tests are passing. This protects you from the day-to-day churn in the
Magpie repository.

### GitLab access

To obtain access an HPC token is required.

## Build Magpie

In a fresh Magpie clone the MyTRIM and SPPARKS submodules need to be checked out
first.

```bash
git submodule update --init
```

Now you are ready to build and test Magpie

```bash
make -j 8
./run_tests -j 8
```

Replace the `8` with the number of CPU cores available on your system. A success
installation should pass all tests and the last lines of output should read.

```text
-------------------------------------------------------------------------------------------------------------
Ran 61 tests in 54.5 seconds.
61 passed, 3 skipped, 0 pending, 0 failed
```
