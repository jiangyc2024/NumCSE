## Numerical Methods for CSE

# Project Structure
<ul>
  <li>The codes used in the script are in <i>LectureCodes</i>.</li>
  <li><i>MathGL</i> contains a documentation and example codes for MathGL,
      as well as the Figure library.</li>
  <li>The CMake Find<Package>.cmake are located in <i>CMake</i></li>
</ul>

# How to use
In the LectureCodes you can find a folder for each code from the script. <br>
There you can find different versions of the code:
<ul>
  <li>a C++ version</li>
  <li>a Matlab (most of the codes)</li>
  <li>a Python (most of the codes)</li>
</ul>
<br>
Building C++-codes & usage with Eigen and MathGL/Figure: <br>
In the LectureCodes directory you can find <code>GetModules.cmake</code>. 
That file defines a function you can use to find the packages 
<code>Eigen</code>, <code>MathGL</code> or <code>Figure</code>.
A typical example of a CMakeLists.txt would be:
<pre><code>
project(myProject)
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



