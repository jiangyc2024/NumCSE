# NumCSE - CMake documentation

## General build procedure

The complete NumCSE project is build from the *CMakeLists.txt* file in the root directory of the project. From here compilation is delegated to all subdirectories that contain code (currently only the `LectureCodes` directory) which then delegate compilation further to their own subdirectories and so on. This is archived either by cmake's built in `add_subdirectory` command or by a convenience macro `add_all_subdirectories` that executes `add_subdirectory` for all subdirectories of the given path (e.g. `CMAKE_CURRENT_SOURCE_DIR`) and warns if no *CMakeLists.txt* file was found in one of them. For project developers this means that in most cases new codes with a *CMakeLists.txt* added to an already existant directory will be build automatically. If they however add new folders with subdirectories, where only the subdirectotires contain code, additional *CMakeLists.txt* files must be added to propagate inclusion of the *CMakeLists.txt* files in those subdirectories using either `add_subdirectory` or `add_all_subdirectories`.

The NumCSE project requires a set of third party libraries that are used in most of the C++ codes. In order to provide a convenient way of writing/including new executables download of the library, if necessary, and setup of header include paths of all third party libraries is done once in the *CMakeLists.txt* of the root directory and linking executables with those libraries is done automatically whenever an executable is created with the `add_executable_numcse` macro.

All compiled binaries are placed in a subdirectory of the `bin` folder (resides in `PROJECT_BINARY_DIR`). The subdirectories are created automatically within the `add_executable_numcse` macro and resamble the path of the calling *CMakeLists.txt* file.

## Macro documentation

### `add_executable_numcse(<name> source1 [source2 ...])`

Adds an executable target to be built from the source files listed in the command invocation, sets the `RUNTIME_OUTPUT_DIRECTORY`, and links with MathGL and the MathGL Figure library. The target name is contrary to cmakes built in `add_executable` command not just the `<name>` parameter, but a combination of the directory path and `<name>` that is determined by the `get_target_name_numcse` macro.

**Example Usage**

	add_executable_numcse(denselinalg main.cpp)

**Note** that due to the different naming of the target created by `add_executable_numcse` linking with libraries must make usage of the `get_target_name_numcse` macro. An example is given below:

	# create executable
	add_executable_numcse(denselinalg main.cpp)

	# link with some_lib
	get_target_name_numcse(denselinalg target_name)
	target_link_libraries(${target_name} some_lib)

### `get_target_name_numcse(<base_name> <output_variable>)`

Assembles a unique target (executable) name from `base_name` and `CMAKE_CURRENT_SOURCE_DIR` and stores the result in `output_variable`.

The target name is assembled by
 1. splitting the current source directory into chunks
 2. where each chunk is converted 
    from CamelCase to hyphen-seperated-lower-case
 3. dropping `Eigen` suffix (if present)
 4. afterwards all chunks are concatinated by underscores.

E.g. calling `get_executable_name_numcse(main)` in `LectureCodes/NumQuad/numquaderrs/Eigen/CMakeLists.txt` will return `lecture-codes_num-quad_main`

**Example Usage**

	# get target name and store it in `target_name`
	get_target_name_numcse(main target_name)
	# under the assumption that the CMakeLists.txt file is located in LectureCodes/NumQuad/
	# the following will print:
	#	lecture-codes_num-quad_main
	message(${target_name})

### `add_all_subdirectories(<directory>)`

Include all *CMakeLists.txt* files in subdirectories of `<directory>` or emit a warning if no such file was found.

**Example Usage**

	add_all_subdirectories(${CMAKE_CURRENT_SOURCE_DIR})