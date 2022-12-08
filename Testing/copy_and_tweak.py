import re
import sys

#
# CURRENT ISSUES:
#


def parseWriteChange(student_sol, copy, tests):
    """
    Checks the c++ file student_sol for function definitions and changes the name of those to *_TEST.
    """
    suffix = "_TEST"
    keyword = "[OUT OF CLASS]"

    # get all function signatures to be tested from the test file
    ncc = False
    function_signatures = []
    for line in tests:
        disable = re.match('// DISABLE_TESTS', line)
        no_class_copy = re.match('// NO_CLASS_COPY', line)
        if no_class_copy:
            ncc = True
        if disable:
            return
        parse = re.match('\s*TEST_CASE\("([\S\s]*)"\s\*.*', line)
        if parse and keyword not in parse.group(1):
            function_signatures.append(parse)
    classes = []

    for file in student_sol:
        # expect include guards, change them
        if file.readline().startswith("#ifndef"):
            next(file)
            next(file)
        elif file.readline().startswith("#pragma once"):
            next(file)
        enum_helper = False
        for line in file:
            write = True
            parse_class = re.match(
                "\s*(class|struct)(\s*)(\S*)(\s*)({|\s*)", line)
            if parse_class:
                copy.write(parse_class.group(1) + parse_class.group(2) + parse_class.group(
                    3) + suffix + parse_class.group(4) + parse_class.group(5) + "\n")
                classes.append(parse_class.group(3))
                write = False
            for el in function_signatures:
                parse = re.match("(\s*)" + "(Eigen::|std::)?" +
                                 re.escape(el.group(1)) + "(\s*\(.*)", line)
                if parse:
                    # parse2 = re.search("(<.*>)$", re.escape(el.group(1)))
                    parse2 = re.search("(<.*>)$", el.group(1))
                    if parse2:
                        copy.write(parse.group(1) + el.group(1).replace(parse2.group(1),
                                   "") + suffix + parse2.group(1) + parse.group(3) + "\n")
                    else:
                        copy.write(parse.group(1) + el.group(1) +
                                   suffix + parse.group(3) + "\n")
                    write = False
            # look for enums and skip the whole enum scope
            if line.startswith("enum") or enum_helper:
                write = False
                enum_helper = True
                if "}" in line:
                    enum_helper = False
                    continue
            # write the new tweaked line but skip endifs and include statements with signature #include "..."
            if write and not line.startswith("#endif") and not line.startswith('#include "'):
                new_line = line
                # replace all class occurences with its test version
                if not ncc:
                    for el in classes:
                        new_line = line.replace(el, el + suffix)
                copy.write(new_line)
        copy.write("\n")


def searchStudentSolution(main):
    """
    Searches the student solution using main.cpp and the exact order of inclusion in it.
    Returns a list of filenames.
    """
    return_list = []
    other_headers = ["polyfit", "polyval", "ode45", "fft", "solution",
                     "solution2", "solution3"]  # other headers ending with .hpp
    for line in main:
        parse = re.match('\s*\#include\s*"(.*).hpp"', line)
        if parse and parse.group(1) not in other_headers:
            return_list.append(parse.group(1) + ".hpp")
    return return_list


def cmake(source_dir, build_dir):
    """
    Will perform the copy and tweak operation if called from a cmake build system inside a build folder.
    The tweaked copy is named 'copy.hpp'.
    Will also put header guards into the file.
    """
    main = open(source_dir + "/template/main.cpp", "r")
    files = searchStudentSolution(main)
    main.close()
    student_sol = []
    for file in files:
        student_sol.append(open(source_dir + "/template/" + file, "r"))

    copy = open(build_dir + "/copy.hpp", "w")
    tests = open(source_dir + "/template/tests.cpp", "r")

    copy.write("#ifndef COPY_HPP\n")
    copy.write("#define COPY_HPP\n")

    for el in files:
        # ../../../../../build/bin/Assignments/PolishedCodes/Introduction/MatrixBlocks/
        copy.write('#include "{}/solution/{}"\n'.format(source_dir, el))

    parseWriteChange(student_sol, copy, tests)

    copy.write("#endif")

    copy.close()
    for item in student_sol:
        item.close()
    tests.close()


def codeexpert():
    """
    Will perform the copy and tweak operation if called from a script that is used by codeexpert.
    The tweaked copy is named 'copy.hpp'.
    Will also put header guards into the file.
    Puts copy.hpp in a writable folder on the server – MUST BE REMOVED AFTER COMPILING!
    """
    main = open("main.cpp", "r")
    files = searchStudentSolution(main)
    main.close()
    student_sol = []
    for file in files:
        student_sol.append(open(file, "r"))

    copy = open("copy.hpp", "w")
    tests = open("tests.cpp", "r")

    copy.write("#ifndef COPY_HPP\n")
    copy.write("#define COPY_HPP\n")

    # complicated but for backward compatibility needed
    copy.write('#include "solution.hpp"\n')
    
    if len(files) > 1:
        for i in range(1, len(files)):
            copy.write('#include "solution{}.hpp"\n'.format(i + 1))

    parseWriteChange(student_sol, copy, tests)

    copy.write("#endif")

    copy.close()
    for item in student_sol:
        item.close()
    tests.close()


if __name__ == "__main__":
    if len(sys.argv) == 3:
        cmake(str(sys.argv[1]), str(sys.argv[2]))
    else:
        codeexpert()
