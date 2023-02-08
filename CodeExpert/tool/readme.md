# Code Expert tool

Use `export.sh` to automate the upload process to CodeExpert. For more information, consult `export.js` and the transcript.

# NCSEFL_Prbrefs.aux

 The file `NCSEFL_Prbrefs.aux` is intentionally tracked in version control, even though it is automatically generated during compilation of LaTeX files in the LaTeX repository for this course. We do this to keep this repository agnostic of any LaTeX files in terms of dependencies. It is best to update `NCSEFL_Prbrefs.aux` when the chapter structure may have changed and assignment codes need to be published to the Code Expert platform. Do this by running `make update_refs` from [npdeflipped/LaTeX_NCSE/Makefile](https://gitlab.math.ethz.ch/ralfh/npdeflipped/-/blob/master/LaTeX_NCSE/Makefile).