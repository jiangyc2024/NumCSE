#Debian GNU/Linux 11 (bullseye)
from rikorose/gcc-cmake
run apt-get update
run apt-get install -y python3 python3-pip python3-dev
run pip3 install numpy

run echo "install clang-tidy 14"
run apt-get install -y lsb-release wget software-properties-common
run wget https://apt.llvm.org/llvm.sh
run chmod +x llvm.sh
run ./llvm.sh 15 all
run ln -s /usr/bin/clang-tidy-15 /usr/bin/clang-tidy

#temporarily cmake on a copy of the repository
copy . ./numcse
run ls
run mkdir /build
workdir /build
run cmake -DTIDY=true /numcse
run make -j7 Eigen
run rm -r /numcse