# Guidance on environment setup using Docker
This Dockerfile could quickly set up an environment that compiles the whole NumCSE project. 

1. Install docker on your computer following https://docs.docker.com/
2. Open a terminal (with administrator authorization or use `sudo`), `cd` into this folder
3. execute `docker build -t numcse .` to build an image
4. execute `docker run -d --cap-add sys_ptrace -p "50:22" --name numcse_env numcse` to make a container based on the image
5. now you can connect into the container using `ssh root@localhost -p 50` with password: `mypassword` (strongly suggest to change that in the Dockerfile)
6. or use any IDE like CLion or VSCode Remote to connect to it
7. you could start or stop the container using `docker start numcse_env` or `docker stop numcse_env`
8. you could delete the container using `docker rm numcse_env` after stopping it and delete the image using `docker rmi numcse`

The use of `ssh` is constrained by CLion since it only support toolchain using ssh. If don't need it, You are free to delete ssh settings in the Docker file and simply use `docker attach` or something else.