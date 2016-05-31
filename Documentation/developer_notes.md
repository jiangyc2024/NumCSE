# Developer notes

This document contains notes relevant for all project developers and can be
ignored by the students.

## Organization of the porting
If you want to port a chapter or a part of it first check if it hasn't 
already been ported. Then see if there exists an 'Issue' for that part.
If not then create a new Issue in which you exactely declare which codes
you will be working on and close you Issue when you're done. <br>

If you're working on a code of the script be sure to copy the Matlab file 
(if it exists) from the SVN repository to the according directory of the code
in the subdirectory 'Matlab', so we keep all in one place. 
If your code is working replace the Matlab code in the script on SVN by your
C++ version.<br>


For now the codes appearing in the script are the priority (and not
codes that only create plots, we're keeping the old ones).

## Editing the `.tex` files
The `.tex` files for the script are organized on SVN.

*Be very careful when editing the Latex code!*

To find the place where to put your code in the script open the 
`NCSENEW_Ch_<yourChapter>.tex` and then use text-search
for the name of the directory of your code, as the directory names
are chosen after the Latex label of the code.

Note: The code environment that's being used is called *samcode*.

## Compiling the `.tex` files

All `.tex` files are included in `NCSENEW.tex`.
So if you changed something in e.g. `NCSENEW_Ch_Evp.tex` you 
can simply compile the whole script to see the changes(*):

	$ latex -interaction=scrollmode -shell-escape NCSENEW.tex

You will be prompted for the 'slidemode', explanations of the modes are given
on the screen just above the question.

In the output `NCSENEW.dvi` the formatting doesn't work very well
so you have to convert it to a `pdf` version with:

	$ dvipdf NCSENEW.dvi NCSENEW.pdf

**(*) you can also compile only your chapter by choosing slidemode 4 and
uncomment the name of the chapter in `NCSENEW.tex` in line 100
and following and comment the other chapters.
However try to compile the whole script when your changes are working on the
small scale to make sure everything still works as it should.**

## Embedding the C++-codes in the script

As the codes are organised on Gitlab and the script on SVN we use a symlink
from `ncse_new/../../`(the parent directory of the SVN-NCSE directory)
to the LectureCodes Gitlab-directory.
To set up the symlink access the parent directory of the SVN-NCSE directory and do:

	ln -s /path/to/Gitlabs/LectureCodes CppScriptCodes

Place the codes in the `.tex` files as follows:

	\begin{lstlisting}
	  \lstinputlisting[style=cpp,&lt;optional arguments&gt;]{%
	    ./../../CppScriptCodes/Chapter/Eigen/code.cpp}
	\end{lstlisting}