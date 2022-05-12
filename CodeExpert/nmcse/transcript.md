## Transcript of the communication with the CodeExpert team

A *[...]* indicates an omitted part.

### 2022/03/17 [Heinrich Grattenthaler](heinrich.grattenthaler@sam.math.ethz.ch):

[...] I am working for Prof. Ralf Hiptmair on the NumCSE course and we are currently rethinking our build process. We have decided that we would like to update both Eigen to 3.4 and the C++ std in our own repository, but naturally, we need to make sure this is compatible with CE in order for students being able to reproduce our codes. Are these updates possible? The CE “expert” in our team told me I would have to contact you (the admins) about this, because Eigen was centrally installed. [...]

###  2022/03/17 [David Sichau](expert@inf.ethz.ch):

[...]

I think you are still using the generic 1 image. As this is based on the outdated redhat-7 we will not make any changes to this image.

We have an generic 2 image which is based on redhat-8 but did not work in some of your use cases.

I think this update might be useful to make the switch to generic 2. We can help you in this process. But the testing of this image must happen on your side.

I gave you access to our environments:
https://gitlab.inf.ethz.ch/OU-LECTURERS/containers/cxenv

There you can have a look at generic-2 image. Please let us know what you actually need to be installed.

For more details on our execution please have a look at our documentation:
https://docs.expert.ethz.ch/lecturers/#runner-and-cxenvironments

[...]

###  2022/03/18 [Heinrich Grattenthaler](heinrich.grattenthaler@sam.math.ethz.ch):

[...]

thank you for your help. I did some investigating about our whole cx workflow, we had an internal talk and it turns out that we have two major issues. I am telling you these, because I believe it is better to give you the big picture.

1. Stable build environments

We have a very large repository of codes maintained by a very heterogeneous team (in terms of machines used in development), thus we are constantly running into build and compatibility issues. The assignment codes that are later uploaded to cx are a part of that repository and they are maintained and developed „offline“, i.e. independent of cx on our machines. To eliminate that risk we are currently working on a unified docker solution to solve all these problems. Since the assignment codes are tested and developed in that environment, it would only make sense that we use this (or some easy-to-maintain mirror of this) container for cx as well. I suggest we develop a custom image / cxEnvironment for our course, (I see this has already been done for some other course called „ifme1“). The drawback is that the NumCSE team will have to maintain the image.

2. Upload automation

Since the exercises have been developed over years and each year only a small subset of assignments are published on cx, we would like to have an automated (in the sense that one cannot make crucial mistakes), scalable way of uploading our codes as a whole project, such that the work necessary in cx itself is limited to some reasonable amount. We have experimented with browser automation (which worked OK for a while), but ultimately this is a costly and unsustainable decision, once again driving up the month-to-month maintenance work for the course. In my opinion, we are in an unfortunate scenario right now. Do you have any ideas on this? Do other courses have similar issues? I am aware that cx was probably not designed with this kind of workflow in mind. Is it a feasible way to maybe convert our file structure to some flat hierarchy which can be manually bulk-uploaded to cx? What if we .zip the directory?

Here is an example of one assignment (which is then mapped to a cx project) and how it is stored in our repository. Our testing scripts are stored externally  and later added to the cx project.

[...]

###  2022/03/19 [Heinrich Grattenthaler](heinrich.grattenthaler@sam.math.ethz.ch):

[...]

I just wanted to apologize for bothering you with the "upload automation“ chapter. Our team must have missed the „upload from file“ feature in cx. This completely solves the second issue in the best way I could ever imagine. Sorry…

[...]

###  2022/03/21 [David Sichau](expert@inf.ethz.ch):

[...]

yes I think a custom images is the way to go for your course as you have very special usecases. Do you want to start with our base image or do you want to base it on generic-2? I will then create a special repo for you. As soon as you are happy with the image let us know and we will deploy it to Code Expert.

The automatic upload was developed based on the needs for your course. I am glad that this solves your problems.

[...]

###  2022/03/25 [Heinrich Grattenthaler](heinrich.grattenthaler@sam.math.ethz.ch):

[...]

I still had to clarify a few issues with our exercise uploading process in the last few days, so I could not yet give you an answer. I think it is better to start fresh with your base image to keep the containers clean. As previously mentioned, we would like to use a single, unified image for Code Expert and internal use cases so that our assistants (and sometimes students) will only have to deal with a single container for all purposes. Do you think this can be achieved or is that a bad idea for some reason I cannot see right now? This means that the image will be larger than it has to be for just the Code Expert use case. An example would be development tools like CMake and Clang-tidy.

[...]

###  2022/03/28 [Stefan Dröschler](expert@inf.ethz.ch):

[...]

we've prepared two new repositories for environment images for you:

[1] https://gitlab.inf.ethz.ch/OU-LECTURERS/containers/cxenv/nmcse
[2] https://gitlab.inf.ethz.ch/OU-LECTURERS/containers/cxenv/nmcse-dev

We suggest you set up and use [1] for the actual student jobs and [2] for development, i.e. [2] can include additional development tools you might need, while [1] is kept as lean as possible. 

Currently only you have been assigned maintainer role in gitlab and pull permissions from our registry cxhub.ethz.ch. Pushing to the repositories will trigger an automatic rebuild of the images (including push to the registry), but will not make them available in Code Expert immediately. Please drop us a mail whenever you have a new version of an image ready, so we can update them in Code Expert accordingly. 

Additionally, if you'd like to change the default project associated with these images, please set one up in Code Expert and send us an export (i.e. "Download Project" in our IDE). 

[...]

###  2022/03/28 [Heinrich Grattenthaler](heinrich.grattenthaler@sam.math.ethz.ch):

[...]

thanks for the setup. I still have to think about the dual repository approach, maybe we will end up basing our dev image on the lean version and just build that part ourselves (but that is still far into the future). To which account/email did you give the pull permissions for cxhub.ethz.ch? I tried two docker accounts registered for heinrich.grattenthaler@sam.math.ethz.ch and heinrich.grattenthaler@inf.ethz.ch, but authentication failed for a pull from cxhub.ethz.ch/cx/cxenv/base-rhel8:latest.

[...]

###  2022/03/28 [Stefan Dröschler](expert@inf.ethz.ch):

[...]

for the moment we think the dual environment approach best suits your use case and our environment. Two repositories are necessary in order to benefit from our automation in the background (we haven't had the need to build multiple environments from one repos so far). You can, of course, base the nmcse-dev image on nmcse. 

As for the permissions: please log in to our registry with your nethz account, i.e. docker login -uhgratten cxhub.ethz.ch

[...]

###  2022/04/03 [Heinrich Grattenthaler](heinrich.grattenthaler@sam.math.ethz.ch):

[...]

I tried building the nmcse image locally. In instruction 4 ``RUN dnf install -y autoconf automake gcc gcc-c++ bison flex b...``, dnf fails to access the isginf-8 package repository: 

```
#8 132.7 isginf Extra Packages for 8 - x86_64            0.0  B/s |   0  B     01:59    
#8 132.7 Errors during downloading metadata for repository 'isginf-8':
#8 132.7   - Curl error (28): Timeout was reached for http://install-el8.inf.ethz.ch/software/el8/x86_64/repodata/repomd.xml [Connection timed out after 30001 milliseconds]
#8 132.7   - Curl error (28): Timeout was reached for http://install-el8.inf.ethz.ch/software/el8/x86_64/repodata/repomd.xml [Connection timed out after 30000 milliseconds]
#8 132.7   - Curl error (28): Timeout was reached for http://install-el8.inf.ethz.ch/software/el8/x86_64/repodata/repomd.xml [Connection timed out after 30002 milliseconds]
#8 132.7 Error: Failed to download metadata for repo 'isginf-8': Cannot download repomd.xml: Cannot download repodata/repomd.xml: All mirrors were tried
```

For now, I ignored the settings in deploy-docker.cfg, but I assume this is unrelated.
Can you help me out here? Maybe I am not supposed to build these images locally at all and only use the CI build for development?

[...]

###  2022/04/04 [Stefan Dröschler](expert@inf.ethz.ch):

[...]

the repository is maintained by ISG and is only accessible internally or via VPN. 

[...]