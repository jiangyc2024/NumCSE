# Links

- Dox:
    https://svn.id.ethz.ch/sam/Numcourses/ncse/ncse_new/documentation/doc_samstyle.pdf
- Macros:
    https://svn.id.ethz.ch/sam/Numcourses/ncse/ncse_new/NCSENEW_macros.tex

# Requirements

We use `unifdef` and `python3` to deploy.

# TODOs

- "lib" mode, see below
- decide variable/functions/class naming convention
- merge tex+cpp?
- ./deploy.cpp must be called within "Deployment", fix this
- TODO: awesome but too time consuming: use executables to automatically
  create binary files and include them in templates/solutions (using json
  description)
- TODO: put `deploy.py`in `make`

# Assignment codes tree structure

The "Assigments" folder is structured in the following way:
 - "./LaTeX": contains all ".tex" files, *one* .tex file per problem. Each
 - ".tex" file is in an appropriate "Chapter" folder.
   Each chapter folder in LaTeX must have a corresponging folder in LectureCodes.
   TODO: possibly merge LaTeX and Codes
 - "/Codes": contains all ".cpp" and other codes that are used for problem sheets.
   Also contains binary data such as generated images.
   Each code goes in an appropriate "Chapter/ProblemName" folder, and must have a corresponding laTex file.
   SPECIAL: each "ProblemName" folder must have a file "description.json"
   containing a list of all files used by the problem. See below.
 - "./Deployment": contains scripts and details for packaging of templates/solutions
 - "./Legacy": svn imported old codes and tex, at some point will be purged

# Writing C++

We use the following conventions:
- try and use one cpp file per problem
- templates and solutions must compile
- give each cpp file a unique name.cpp
- for range based listing: use:

    ```
    /* SAM_LISTING_BEGIN_3 */
    // CODE Goes here
    /* SAM_LISTING_END_3 */
    ```

    WARNING: due to issues with lstlisting, we have the following restrictions:
    - code between multiline comments is not escaped, must compile with latex
    - TODO: escape code

# Mark solution-only code within:

    ```
    #if SOLUTION
    // CODE Goes here
    #endif \\ SOLUTION
    ```

# Mark internal-only code within:

    ```
    #if INTERNAL
    // CODE Goes here
    #endif \\ INTERNAL
    ```

# Include code into LaTeX files:

Use the range based lstinput listing. To include template of file "file.cpp"
include file as if it were "templates/file.cpp" instead. For the
solution, include the file as "solutions/file.cpp" instead.

`deploy.py` will generate those files for you
TODO: put deploy in `make`

# Create description.json

Each problem directory has a "description.json" file. This contains 5 fields:
- "name" a unique name for the problem
- "templates" a list of files included as templates into the bundle
- "solutions" a list of files included as solutions into the bundle
- "all" a list of files included everywhere
- "internal" a list of files not bundled (used for LaTeX and inside repository only)
- TODO: "lib" a library similar to "include" in the main JSON,
  but will be copied only o the selected problem 

# Update the assigment bundle list (edit assignment_list.json)

1. Add a new problem sheet by adding a new entry into the "ProblemSheet" array.
2. Give a name for the problem sheet and a list of files in "problems".
3. The path is relative to the "assignment_dir" path.
4. The "working dir" is where all bundles will be created.
5. The "cmake" entry specifies the cmake header to use (do not edit).
6. The "include" field specifies files and folders that are included in every bundle.
7. "Orphan" are problems not (yet?) included in a problem sheet.

# Deploy

Call:

```
./deploy.py
```

To deploy newly created C++ files. Within each "problem" folder, 4 new folders
are created:
- "templates"
- "solutions"
- "solutions_nolabels"
- "tempaltes_nolabels"

Additionally in `working_dir` a folder is created for each "Problem Sheet".
Folder are also compressed. These folders contain solutions and templates
with labels stripped away.

# Conventions:

Follow the following rules:
- use least amount of templates
- we use `using Eigen`, we do not `using std`
- function name: TODO (maybe `underscore_separated_name`)
