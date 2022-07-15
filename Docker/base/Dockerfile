from gcc:11

# ==== install cmake ===============================================================

# based on https://github.com/Rikorose/gcc-cmake/blob/1e2943e709239beed1d368bf5c511823f956b54e/Dockerfile
arg cmake_version=3.23.2

run platform="$(uname -m)" \
    && wget https://github.com/Kitware/CMake/releases/download/v${cmake_version}/cmake-${cmake_version}-linux-$platform.sh \
        -q -O /tmp/cmake-install.sh \
	&& chmod u+x /tmp/cmake-install.sh \
    && mkdir /usr/bin/cmake \
    && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
    && rm /tmp/cmake-install.sh

env PATH="/usr/bin/cmake/bin:${PATH}"

# ==== install python ==============================================================

run apt-get update
run apt-get install -y python3 python3-dev python3-numpy python3-matplotlib

# ==== install clang-tidy ==========================================================

arg clang_tidy_version=15

run apt-get update
run apt-get install -y lsb-release wget software-properties-common
run wget https://apt.llvm.org/llvm.sh \
    && chmod +x llvm.sh \
    && ./llvm.sh ${clang_tidy_version} all \
    && ln -s /usr/bin/clang-tidy-${clang_tidy_version} /usr/bin/clang-tidy

# ==== install Eigen ===============================================================

arg eigen_version=3.4.0 

run apt-get update
run apt-get -y install unzip
run mkdir -p "/eigen" cd "/eigen" \
      && curl -L https://gitlab.com/libeigen/eigen/-/archive/${eigen_version}/eigen-${eigen_version}.zip > ./eigen.zip \
      && unzip eigen.zip \
      && mv eigen-${eigen_version} /usr/include/eigen3 \
      && rm -rf "/eigen" \
      || exit 1