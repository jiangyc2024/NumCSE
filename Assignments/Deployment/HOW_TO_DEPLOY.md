* Requirements

unifdef and python3

* Assignment codes tree structure

The "Assigments" folder is structured in the following way:
 - "/LaTeX": contains all ".tex" files, *one* .tex file per problem. Each ".tex" file is in an appropriate "Chapter/ProblemName" folder.
   Each folder in LaTeX must have a corresponging folder in LectureCodes.
   TODO: possibly merge LaTeX and Codes
 - "/Codes": contains all ".cpp" and other codes that are used for problem sheets. Also contains binary data such as generated images.
   Each code goes in an appropriate "Chapter/ProblemName" folder, and must have a corresponding laTex file.
   SPECIAL: each ProblemName folder must have a file "description.json" containing a list of all files used by the problem. See below.
 - "./Deployment": contains scripts and details for packaging of tempaltes/solutions

 - "Legacy": svn imported old codes and tex

* Writing C++

We use the following conventions:
- try and use one cpp file per problem
- give each cpp file a unique name.cpp
- for range based listing: use:

```
/* SAM_LISTING_BEGIN_3 */
// CODE Goes here
/* SAM_LISTING_END_3 */
```

- mark solution-only code within:

```
#if SOLUTION
// CODE Goes here
#endif \\ SOLUTION
```

- mark internal-only code within:

```
#if INTERNAL
// CODE Goes here
#endif \\ INTERNAL
```

* Include code into LaTeX files:

use the range based lstinput listing. To include template of file "file.cpp" include file "templates/file.cpp" instead. For
solution include file "solutions/file.cpp" instead.

* Create description.json

Each problem directory has a "description.json" file. This contains 5 fields:
- "name" a unique name for the problem
- "templates" a list of files included as templates into the bundle
- "solutions" a list of files included as solutions into the bundle
- "all" a list of files included everywhere
- "internal" a list of files not bundled (used for LaTeX and inside repository only)

* Update the assigment bundle list (edit assignment_list.json)

Add a new problem sheet by adding a new entry into the "ProblemSheet" array. Give a name for the problem sheet and a list of
files in "problems". The path is relative
to the "assignment_dir" path. The working dir is where all bundles will be created. The "cmake" entry specifies the cmake header to use (do not edit).
The "include" field specifies files and folders that are included in every bundle.

* Deploy

Call:

```
./deploy.py
```
