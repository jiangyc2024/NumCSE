#https://blog.jetbrains.com/clion/2020/01/using-docker-with-clion/
#build: docker build -t numcse .
#start: docker run -d --cap-add sys_ptrace -p "50:22" --name numcse_env numcse
#NOTE: port 50 on localhost is opened for ssh access to the container
FROM ubuntu:18.04

RUN apt-get update
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt-get install -y git cmake libpng++-dev freeglut3-dev libboost-all-dev zlib1g zlib1g-dev \
    openssh-server build-essential gcc g++ gdb clang rsync tar python3 python3-dev python3-numpy python3-matplotlib \
    libfftw3-3 libfftw3-dev \
  && apt-get clean

RUN mkdir /var/run/sshd
RUN echo 'root:mypassword' | chpasswd
RUN sed -i 's/.*PermitRootLogin.*/PermitRootLogin yes/g' /etc/ssh/sshd_config \
  && sed -i 's/.*PasswordAuthentication.*/PasswordAuthentication yes/g' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

RUN git clone https://gitlab.math.ethz.ch/NumCSE/NumCSE.git

EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]

RUN wget https://github.com/Kitware/CMake/releases/download/v3.17.0/cmake-3.17.0-Linux-x86_64.sh \
  && mkdir /opt/cmake \
  && sh /cmake-3.17.0-Linux-x86_64.sh --prefix=/opt/cmake --skip-license \
  && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

