FROM gcc:11

# ==== install cmake ===============================================================

# based on https://github.com/Rikorose/gcc-cmake/blob/1e2943e709239beed1d368bf5c511823f956b54e/Dockerfile
ARG CMAKE_VERSION=3.23.2

RUN PLATFORM="$(uname -m)" \
	  && wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-$PLATFORM.sh \
	  -q -O /tmp/cmake-install.sh \
	  && chmod u+x /tmp/cmake-install.sh \
      && mkdir /usr/bin/cmake \
      && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
      && rm /tmp/cmake-install.sh

env PATH="/usr/bin/cmake/bin:${PATH}"

# ==== install python =================================================================

ARG PYTHON_VERSION=3.7.3

RUN apt-get update
RUN apt-get install -y build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev curl libbz2-dev
RUN mkdir -p "/python" && cd "/python" \
 	  && curl https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz > ./Python-${PYTHON_VERSION}.tgz \
 	  && tar -xvf ./Python-${PYTHON_VERSION}.tgz \
 	  && ./Python-${PYTHON_VERSION}/configure --enable-optimizations \
 	  && make altinstall \
      && cd "/" && rm -rf "/python" \
      # create symlinks: python, python3, pip and pip3
      && ln -s /usr/local/bin/python${PYTHON_VERSION%.*} /usr/local/bin/python3 \
      && ln -s /usr/local/bin/python${PYTHON_VERSION%.*} /usr/local/bin/python \
      && ln -s /usr/local/bin/pip${PYTHON_VERSION%.*} /usr/local/bin/pip3 \
      && ln -s /usr/local/bin/pip${PYTHON_VERSION%.*} /usr/local/bin/pip \
      # install pip and all packages
      && pip install --upgrade pip \
      && pip install matplotlib \
      || exit 1

# ==== install clang-tidy ==========================================================

run apt-get install -y lsb-release wget software-properties-common
run wget https://apt.llvm.org/llvm.sh
run chmod +x llvm.sh
run ./llvm.sh 15 all
run ln -s /usr/bin/clang-tidy-15 /usr/bin/clang-tidy