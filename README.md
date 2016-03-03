## Numerical Methods for CSE

# Project Structure
<ul>
  <li>The codes used in the script are in <i>LectureCodes</i>.</li>
  <li><i>MathGL</i> contains a documentation and example codes for MathGL,
      as well as the Figure library.</li>
  <li>The Find&lt;Package&gt;.cmake are located in <i>CMake</i></li>
</ul>

### Organization of the porting
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

### Editing the <code>.tex</code> files
The <code>.tex</code> files for the script are organized on SVN. <br><br>
<i>Be very careful when editing the Latex code!</i> <br><br>
To find the place where to put your code in the script open the 
<code>NCSENEW_Ch_&lt;yourChapter&gt;.tex</code> and then use text-search
for the name of the directory of your code, as the directory names
are chosen after the Latex label of the code. <br>
Note: The code environment that's being used is called <i>samcode</i>.

### Compiling the <code>.tex</code> files
All <code>.tex</code> files are included in <code>NCSENEW.tex</code>.
So if you changed something in e.g. <code>NCSENEW_Ch_Evp.tex</code> you 
can simply compile the whole script to see the changes(*):
<pre><code>$ latex -interaction=scrollmode -shell-escape NCSENEW.tex</code></pre>
You will be prompted for the 'slidemode', explanations of the modes are given
on the screen just above the question. <br>
In the output <code>NCSENEW.dvi</code> the formatting doesn't work very well
so you have to convert it to a <code>pdf</code> version with:
<pre><code>$ dvipdf NCSENEW.dvi NCSENEW.pdf</code></pre> 

<h6>(*) you can also compile only your chapter by choosing slidemode 4 and
uncomment the name of the chapter in <code>NCSENEW.tex</code> in line 100
and following and comment the other chapters.
However try to compile the whole script when your changes are working on the
small scale to make sure everything still works as it should.</h6>

# How to use
In the LectureCodes you can find a folder for each code from the script. <br>
There you can find different versions of the same code:
<ul>
  <li>a C++ version</li>
  <li>a Matlab (most of the codes)</li>
  <li>a Python (most of the codes)</li>
</ul>

### Executing C++-codes:
Every C++-code that is supposed to be executed (i.e. to create plots) should
contain a CMakeLists.txt.
Then you can do:
<pre><code>$ mkdir build && cd build 
$ cmake ..
$ make 
$ ./exec_name 
</code></pre>


### Building C++-codes & usage with Eigen and MathGL/Figure: <br>
In the LectureCodes directory you can find <code>GetModules.cmake</code>. 
That file defines a function you can use to find the packages 
<code>Eigen</code>, <code>MathGL</code> or <code>Figure</code>.
A typical example of a CMakeLists.txt would be:
<pre><code>project(myProject)
cmake_minimum_required(VERSION 2.8)

# be sure to include the correct path!
include(../GetModules.cmake) 

# this will include the Eigen and Figure directories in DIRS
# and the libraries in LIBS
get_modules("Eigen Figure") 
include_directories(${DIRS})
add_executable(main main.cpp)
target_link_libraries(${LIBS})
</code></pre>


