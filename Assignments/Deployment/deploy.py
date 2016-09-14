#! /usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Filippo Leonardi"
__email__ = "filippo.leonardi@sam.math.ethz.ch"

import os
import json
import subprocess
import tempfile
import zipfile
import shutil

# CMake string to use for each problem
cmake_problem_template = "\
# Make rules for {problem_name}\
"

# String to use for each C++ file in the CMakeLists.txt
# cpp_name_noext is the name of the exeutable
# cpp_name is the name of the .cpp file
cmake_template = "\
# make {cpp_name_noext}\n\
set(SRCS {problem_name}/{cpp_name})\n\
add_executable({cpp_name_noext} ${{SRCS}})\n\n\
"

# String to use for each C++ file in the CMakeLists.txt
# cpp_name_noext is the name of the exeutable
# cpp_name is the name of the .cpp file
cmake_single_problem_template = "\
# make {cpp_name_noext}\n\
set(SRCS {cpp_name})\n\
add_executable({cpp_name_noext} ${{SRCS}})\n\n\
"

def generate_cmake(list_of_cpp_files,
                   outdir, header_file, problem_sheet):
    """
    Given a list of cpp files and an output directory, write a CMakeLists.txt file
    using the header provided in header file. problem_shhet is the name of the
    problem sheet.
    """
    print("Generating CMake file for '{}'.".format(problem_sheet))
    infile = open(header_file, "r")

    outfile = open(outdir + "CMakeLists.txt", "w")

    for line in infile.readlines():
        outfile.write(line.format(problem_sheet_name=problem_sheet))

    for problem in list_of_cpp_files:
        outfile.write(cmake_problem_template.format(problem_name=problem))
        for cpp_file in list_of_cpp_files[problem]:
            cpp_file_noext = ".".join(cpp_file.split(".")[0:-1])
            outfile.write(cmake_template.format(problem_name=problem,
                                                cpp_name_noext=cpp_file_noext, cpp_name=cpp_file))

    infile.close()
    outfile.close()
    
    
def generate_cmake_for_problem(list_of_cpp_files,
                   outdir, header_file, problem_name):
    """
    Given a list of cpp files and an output directory, write a CMakeLists.txt file
    using the header provided in header file. problem_shhet is the name of the
    problem sheet.
    """
    print("Generating CMake file for '{}'.".format(problem_name))
    infile = open(header_file, "r")

    outfile = open(outdir + "CMakeLists.txt", "w")

    for line in infile.readlines():
        outfile.write(line.format(problem_sheet_name=problem_name))

    for cpp_file in list_of_cpp_files:
        cpp_file_noext = ".".join(cpp_file.split(".")[0:-1])
        outfile.write(cmake_single_problem_template.format(problem_name=problem_name,
                     cpp_name_noext=cpp_file_noext, cpp_name=cpp_file))

    infile.close()
    outfile.close()

def copy(file, indir, outdir):
    """
    Aggressively copy file to outdir. File can be a directory and can
    already exists.
    """
    print("Moving '{}' to '{}'".format(indir + file, outdir + file))
    if file == "":
        shutil.copytree(indir, outdir)
    elif file[-1] == "/":
        try:
            shutil.rmtree(outdir + file)
        except:
            pass
        shutil.copytree(indir + file, outdir + file)
    else:
        shutil.copy2(indir + file, outdir + file)

def deploy(filename, indir, outdir, header = None,
           with_solution = False, with_internal = False):
    """
    Parse a C++ file trough "unifdef" and remove SOLUTIONS and INTERNAL ifdef
    code blocks. Must have "unifdef" program. Can toggle if INTERNAL and SOLUTIONS are
    true or false.
    """
    if is_cpp(indir + filename):
        cmd = ["unifdef"]
        cmd.append(indir + filename)
        if with_internal:
            cmd.append("-DINTERNAL=1")
        else:
            cmd.append("-DINTERNAL=0")
        if with_solution:
            cmd.append("-DSOLUTION=1")
        else:
            cmd.append("-DSOLUTION=0")

        cmd.append("-o")
        cmd.append(outdir + filename)

        print("Preparing solutions/templates for '{}'".format(filename))
        print(" ".join(cmd))
        subprocess.call(cmd)
        
        if header:
            f = open(outdir + filename, "r")
            lines = []
            for line in f.readlines():
                lines.append(line)
            f.close()
            
            f = open(outdir + filename, "w")
            f.write(header)
            for line in lines:
                f.write(line)
            f.close()
    else:
        copy(filename, indir, outdir)

def strip_labels(filename, indir, outdir = None):
    """
    Remove labels from files. Currently rmoves "SAM_LISTING" comment blocks.
    """

    file = open(indir + filename, "r")

    if outdir == None:
        split_file = filename.split(".")
        split_file[-2] += "_strip"
        out_filename = ".".join(split_file)
    else:
        out_filename = outdir + filename.split("/")[-1]

    print("Stripping '{}', writing to '{}'.".format(filename, out_filename))
    out_file = open(out_filename, "w")

    for line in file.readlines():
        if "SAM_LISTING" not in line:
            out_file.write(line)

    file.close()
    out_file.close()

def mkdir(dir):
    """
    Aggressively create a directory if non existing.
    """
    if not os.path.exists(dir):
        os.makedirs(dir)

def is_cpp(file):
    """
    Returns True if file looks like a C++ file (header of .cpp)
    """
    return "cpp" == file.split(".")[-1]

def open_problem_description(assignment_dir, problem_obj, main_obj):

    print("---> Processing '{}'".format(problem_obj))

    problem_dir = assignment_dir + problem_obj

    assert(os.path.isdir(problem_dir))

    f = open(problem_dir + "description.json", 'r')

    try:
        obj = json.load(f)
    except json.decoder.JSONDecodeError as err:
        print("Malformed JSON '{}'".format(problem_dir + "description.json"))
        print(err)
        exit(-1)

    if not problem_dir.endswith("/"):
        problem_dir += "/"
        
    try: author = problem_obj["author"]
    except: author = main_obj["author"]
    try: contributors = problem_obj["contributors"]
    except: contributors = main_obj["contributors"]
    try: email = problem_obj["email"]
    except: email = main_obj["contributors"]

    header = open(main_obj["working_dir"] + main_obj["cpp_header"], 'r')
    header_string = ""
    for line in header.readlines():
        header_string += line.format(author=author, 
                                     contributors=contributors, 
                                     email=email)
    header.close()

    f.close()

    return [problem_dir, obj, header_string]

def generate_templates_and_solutions(assignment_dir, problem_dir, obj, header):
    mkdir(problem_dir + "solutions/")
    mkdir(problem_dir + "templates/")

    shared_list = []
    try:
        shared_list += obj["all"]
    except: pass
    try:
        shared_list += obj["internal"]
    except: pass

    if "solutions" in obj:
        for solution_file_path in obj["solutions"]:
            deploy(solution_file_path, problem_dir, problem_dir + "solutions/", header, True)
    if "templates" in obj:
        for template_file_path in obj["templates"]:
            deploy(template_file_path, problem_dir, problem_dir + "templates/", header, False)
    for shared_file_path in shared_list:
        deploy(shared_file_path, problem_dir, problem_dir + "solutions/", header, True)
        deploy(shared_file_path, problem_dir, problem_dir + "templates/", header, False)

def generate_nolabels(assignment_dir, problem_dir, obj, main_obj):

    mkdir(problem_dir + "solutions_nolabels/")
    mkdir(problem_dir + "templates_nolabels/")

    git_add_template_list = []
    git_add_solution_list = []
    template_cpp_list = []
    solution_cpp_list = []

    if "solutions" in obj:
        for solution_file_path in obj["solutions"]:
            file = problem_dir + "solutions/" + solution_file_path
            if is_cpp(file):
                strip_labels(solution_file_path, problem_dir + "solutions/", problem_dir + "solutions_nolabels/")
                solution_cpp_list.append(solution_file_path)
            else:
                copy(solution_file_path, problem_dir + "solutions/", problem_dir + "solutions_nolabels/")
            git_add_solution_list.append(problem_dir + "solutions_nolabels/" + solution_file_path)
    if "templates" in obj:
        for template_file_path in obj["templates"]:
            file = problem_dir + "templates/" + template_file_path
            if is_cpp(file):
                strip_labels(template_file_path, problem_dir + "templates/", problem_dir + "templates_nolabels/")
                template_cpp_list.append(template_file_path)
            else:
                copy(template_file_path, problem_dir + "templates/", problem_dir + "templates_nolabels/")
            git_add_template_list.append(problem_dir + "templates_nolabels/" + template_file_path)
    if "all" in obj:
        for shared_file_path in obj["all"]:
            file = problem_dir + "solutions/" + shared_file_path
            if is_cpp(file):
                strip_labels(shared_file_path, problem_dir + "solutions/", problem_dir + "solutions_nolabels/")
                strip_labels(shared_file_path, problem_dir + "templates/", problem_dir + "templates_nolabels/")
                template_cpp_list.append(shared_file_path)
                solution_cpp_list.append(shared_file_path)
            else:
                copy(shared_file_path, problem_dir + "templates/", problem_dir + "templates_nolabels/")
                copy(shared_file_path, problem_dir + "solutions/", problem_dir + "solutions_nolabels/")
            git_add_template_list.append(problem_dir + "templates_nolabels/" + shared_file_path)
            git_add_solution_list.append(problem_dir + "solutions_nolabels/" + shared_file_path)
                
    generate_cmake_for_problem(template_cpp_list, problem_dir + "templates_nolabels/", main_obj["cmake"] , obj["name"])
    generate_cmake_for_problem(solution_cpp_list, problem_dir + "solutions_nolabels/", main_obj["cmake"] , obj["name"])
    git_add_solution_list.append(problem_dir + "solutions_nolabels/CMakeLists.txt")
    git_add_template_list.append(problem_dir + "templates_nolabels/CMakeLists.txt")
    for file in main_obj["include"]:
        copy(file, main_obj["working_dir"], problem_dir + "templates_nolabels/")
        copy(file, main_obj["working_dir"], problem_dir + "solutions_nolabels/")
        git_add_solution_list.append(problem_dir + "solutions_nolabels/" + file)
        git_add_template_list.append(problem_dir + "templates_nolabels/" + file)
    
    subprocess.call(["git", "add"] + git_add_template_list)
    subprocess.call(["git", "add"] + git_add_solution_list)
    #os.subprocess()

    return [None, template_cpp_list, None, solution_cpp_list]

def parse_json(filename):
    """
    Parse a JSON destription of the Assignment bundles.
    """

    solutions_folder_name = "solutions/"
    templates_folder_name = "templates/"
    solutions_nolabels_folder_name = "solutions_nolabels/"
    templates_nolabels_folder_name = "templates_nolabels/"

    archive_format = "zip"
    
    this_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(this_dir)

    f = open(filename, 'r')

    obj = json.load(f)
    

    assignment_dir = obj["assignment_dir"]
    working_dir = obj["working_dir"]
    if not assignment_dir.endswith("/"):
        assignment_dir += "/"
    if not working_dir.endswith("/"):
        working_dir += "/"

    print("Looking for files in '{}'".format(assignment_dir))
    print("Writing output to: {}".format(working_dir))

    for problem_sheet_obj in obj["ProblemSheets"]:
        problem_sheet = problem_sheet_obj["name"]

        problem_sheet_solutions_dir = working_dir + problem_sheet + "/" + solutions_folder_name
        problem_sheet_templates_dir = working_dir + problem_sheet + "/" + templates_folder_name

        try: shutil.rmtree(working_dir + problem_sheet)
        except: pass

        template_cpp_list = {}
        solution_cpp_list = {}

        for problem in problem_sheet_obj["problems"]:
            [problem_dir, problem_obj, header_str] = open_problem_description(assignment_dir, problem, obj)

            try: shutil.rmtree(problem_dir + solutions_nolabels_folder_name)
            except: pass
            try: shutil.rmtree(problem_dir + templates_nolabels_folder_name)
            except: pass

            generate_templates_and_solutions(assignment_dir, problem_dir, problem_obj, header_str)

            [template_file_list, problem_template_cpp_list,
             solution_file_list, problem_solution_cpp_list] = generate_nolabels(assignment_dir, problem_dir, problem_obj, obj)

            copy("", problem_dir + solutions_nolabels_folder_name, problem_sheet_solutions_dir + problem_obj["name"])
            copy("", problem_dir + templates_nolabels_folder_name, problem_sheet_templates_dir + problem_obj["name"])

            template_cpp_list[problem_obj["name"]] = problem_template_cpp_list
            solution_cpp_list[problem_obj["name"]] = problem_solution_cpp_list

        print("\
===============\n\
Bundling {}\n\
===============\
".format(problem_sheet))

        if template_cpp_list != []:
            generate_cmake(template_cpp_list, problem_sheet_templates_dir, obj["cmake"], problem_sheet)
        if solution_cpp_list != []:
            generate_cmake(solution_cpp_list, problem_sheet_solutions_dir, obj["cmake"], problem_sheet)

        for file in obj["include"]:
            copy(file, working_dir, problem_sheet_solutions_dir)
            copy(file, working_dir, problem_sheet_templates_dir)

        print("Creating archive for '{}'".format(problem_sheet))
        shutil.make_archive(problem_sheet + "_template", archive_format, problem_sheet_templates_dir)
        shutil.make_archive(problem_sheet + "_solution", archive_format, problem_sheet_solutions_dir)

    for problem in obj["Orphan"]:

        [problem_dir, problem_obj, header_str] = open_problem_description(assignment_dir, problem, obj)

        try: shutil.rmtree(problem_dir + solutions_nolabels_folder_name)
        except: pass
        try: shutil.rmtree(problem_dir + templates_nolabels_folder_name)
        except: pass

        generate_templates_and_solutions(assignment_dir, problem_dir, problem_obj, header_str)

if __name__ == "__main__":
    
    try: subprocess.call(["unifdef", "--help"])
    except FileNotFoundError:
        print("You must install 'unifdef'.")
        exit(-1)
    
    parse_json("assignment_list.json")
