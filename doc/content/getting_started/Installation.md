# Installation

## Install MOOSE

Refer to the [MOOSE Getting Started](http://mooseframework.org/getting-started/)
pages to set up MOOSE System Environment. Do not clone MOOSE yet, it will be
distributed as a Magpie git submodule.

## Install Magpie

Magpie is hosted on [GitHub](https://github.com/idaholab/magpie)
and can be cloned directly from there using [Git](https://git-scm.com/). We recommend
creating a directory named `projects` to put all of your MOOSE related work (such
as Magpie) in which leads to the following commands (from your home directory):

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

In a fresh Magpie clone the MOOSE and libMesh submodules need to be checked out
and built first.

```bash
git submodule update --init
cd moose
./scripts/update_and_rebuild_libmesh.sh
cd ..
```

Now you are ready to build and test Magpie

```bash
make -j 8
./run_tests -j 8
```

Replace the `8` with the number of CPU cores available on your system.
A success installation should pass all tests and the last lines of output should
read.

```text
-------------------------------------------------------------------------------------------------------------
Ran 48 tests in 24.5 seconds
48 passed, 3 skipped, 0 pending, 0 failed
```
