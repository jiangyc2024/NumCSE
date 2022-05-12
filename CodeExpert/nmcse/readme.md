# Code Expert & nmcse

This directory is meant for development of our custom *nmcse* image on Code Expert. The actual *nmcse* repository is hosted by *isg.inf.ethz.ch* [here](https://gitlab.inf.ethz.ch/OU-LECTURERS/containers/cxenv/nmcse). Use this for publishing changes.

## Basics

Code Expert (cx) is used for our student assignments. It uses Docker containers to run each student's code in an isolated environment. [See documentation](https://docs.expert.ethz.ch/lecturers/).

## History

We used to rely on the cx's own *generic-1* environment/image, but in February 2022 we decided to upgrade our Eigen version to 3.4 and this meant we had to start maintaining our own image. Generally, it is a good idea, since we have quite the unique requirements and it is best we define those ourselves. For more information on how everything came to be, consult the communication transcript.

## Deployment

Pushing to the repository mentioned above triggers an automatic cloud-rebuild of the image . It is then available via the container registry `cxhub.ethz.ch/cx/cxenv/nmcse`. However, to actually make it available on Code Expert, you have to contact expert@inf.ethz.ch. Tell them that you have a new version ready and ask them to update in Code Expert.

## Running

Testing whether the image fulfills our requirements is a bit tricky. In general, you will want to manually pick a few assignments and make sure they can be compiled, run, etc. There are 3 potential ways to test:

* Locally, using the locally built image.
* Locally, using the cloud-built registry image (should not be much different)
* On the cx platform itself, selecting the image when you create a new task.

The latter is the most accurate/true-to-real-conditions test, since it is exactly what the students will use, but it has a very slow iteration speed, obviously.

### Local Build

For building and the image yourself, you have to `docker login -u<nethz> cxhub.ethz.ch` to the code expert registry, which enables you to pull the base image. You may have to ask the CodeExpert team for access to this registry first. Also, you need to be inside the ETH VPN, because the build process requires access an internal package repository.

`cd` to your local clone of the *nmcse* repository and run the `build.sh` next to **this** file.

### Cloud Build

If you want to use the cloud-built image, (after pushing to the *nmcse* repository) run `pull.sh`.

### Local testing

Using either one of the two methods above, you will now have a local image tagged `cx-nmcse`. Test the image by running `run.sh`. This will run an action (either "compile", "run", or "test") on the project in `projectfiles`. You can change the action in the script. To use a different assignment, a directory `projectfiles` (next to **this** file) must be supplied with the assignment data in the same structure that Code Expert projects use. For this, you can export a known-to-work task from the Code Expert platform and pick one of the projects (either "solution" or "template"). Then paste those into `projectfiles`.

**Warning**: Since this directory is bind mounted, changes to this directory will persist on the host machine. This is good for debugging purposes, but may corrupt the files if something goes wrong. It is recommended to keep a clean backup of the `projectfiles` somewhere.

### Cleaning up

Once in a while, you will want to clean up excess images and containers, especially during container development. For this, use `clean.sh`. This is a tentative clean-up process. For a more agressive version, consult `docker image ls` and `docker image rm`.

## Scripts directory

As of now, we do not make use of the default `scripts` directory for Code Expert environments, which can be specified for an image. We leave those to the default. The actions in the `scripts` directory will be overridden by our assignments.
