#
# CODE EXPERT : NUMCSE CX_ENVIRONMENT
#
# This Dockerfile defines your Code Expert cxEnvironment.
# It inherits the settings of the cxenv/base image and
# defines your custom setup.
#
# On top of that, create scripts for the actions `compile`,
# `run`, `test` and `submit` (only the one your cxEnvironment
# supports) and place them in a directory named `script` within
# your Environment's directory.


# Select one of the cxrun base images.
#from cxhub.ethz.ch/cx/cxenv/base-rhel8:latest
from registry.access.redhat.com/ubi8/ubi:8.1

# ==== install essentials ==========================================================

run dnf -y install clang make

# ==== install python ==============================================================

run dnf -y install python3-devel \
      && pip3 install --upgrade pip \
      && pip3 install --upgrade Pillow \
      && pip3 install matplotlib

# ==== install cmake ===============================================================

run dnf -y install wget

# based on https://github.com/Rikorose/gcc-cmake/blob/1e2943e709239beed1d368bf5c511823f956b54e/Dockerfile
arg CMAKE_VERSION=3.22.3

run wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh \
      && chmod u+x /tmp/cmake-install.sh \
      && mkdir /usr/bin/cmake \
      && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
      && rm /tmp/cmake-install.sh

env PATH="/usr/bin/cmake/bin:${PATH}"

# ==== install clang-tidy ==========================================================

run dnf -y install xz
# pre-built binary clang-tidy-15.0.0-x86_64-linux-rhel-8.1.tar.xz 
run wget -O "clang-tidy.tar.xz" "https://polybox.ethz.ch/index.php/s/YSUkzDatd8hdrh4/download" \
      && tar -xf "clang-tidy.tar.xz" -C /usr/local

env PATH="/usr/local/llvm/bin:${PATH}"

# ==== prepare repository ==========================================================

# Temporarily cmakes on a copy of the repository and compile Eigen before everything 
# else. Since parallel compilation fails if Eigen is not installed, this is the best
# way to compile the entire repository as fast as possible.

copy . ./numcse
arg ncpus
run mkdir /build \
      && cd /build \
      && cmake -DTIDY=1 /numcse \
      && make -j${ncpus} Eigen \
      && rm -r /numcse

workdir /build

# TODO: use the actual v15 headers here (llvm-project/clang/lib/Headers)
#copy CodeExpertTool/Headers /usr/lib64/clang/15.0.0/include

# Set the labels for your cxEnvironment.
LABEL cxrun.name="NumCSE (Version 1.0)" \
      maintainer="NumCSE coding team <jinghao.cao@sam.math.ethz.ch>" \
      description="Code Expert NumCSE cxEnvironment"

# Set the USER and WORKDIR directive to its according ENV
# variable. This needs to be at the end of the file.
USER ${USER}
WORKDIR ${WORKDIR}
