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
	function_signatures = []
	for line in tests:
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
		for line in file:
			write = True
			parse_class = re.match("\s*(class|struct)(\s*)(\S*)(\s*)({|\s*)", line)
			if parse_class:
				copy.write(parse_class.group(1) + parse_class.group(2) + parse_class.group(3) + suffix + parse_class.group(4) + parse_class.group(5) + "\n")
				classes.append(parse_class.group(3))
				write = False
			for el in function_signatures:
				parse = re.match("(\s*)" + "(Eigen::|std::)?" + re.escape(el.group(1)) + "(\s*\(.*)", line)
				if parse:
					# parse2 = re.search("(<.*>)$", re.escape(el.group(1)))
					parse2 = re.search("(<.*>)$", el.group(1))
					if parse2:
						copy.write(parse.group(1) + el.group(1).replace(parse2.group(1), "") + suffix + parse2.group(1) + parse.group(3) + "\n")
					else:
						copy.write(parse.group(1) + el.group(1) + suffix + parse.group(3) + "\n")
					write = False
			# if write and not line.startswith("#endif") and not line.startswith("#include") and not line.startswith("enum"): # skip all one liner enums
			if write and not line.startswith("#endif") and not line.startswith('#include "') and not line.startswith("enum"): # skip all one liner enums
				new_line = line
				# replace all class occurences with its test version
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
	other_headers = ["polyfit", "polyval", "ode45", "fft", "solution", "solution2", "solution3"] # other headers ending with .hpp
	for line in main:
		parse = re.match('\s*\#include\s*"(.*).hpp"', line)
		if parse and parse.group(1) not in other_headers:
			return_list.append(parse.group(1) + ".hpp")
	return return_list

def cmake(path_string):
	"""
	Will perform the copy and tweak operation if called from a cmake build system inside a build folder.
	The tweaked copy is named 'copy.hpp'.
	Will also put header guards into the file.
	"""
	main = open(path_string + "/main.cpp", "r")
	files = searchStudentSolution(main)
	main.close()
	student_sol = []
	for file in files:
		student_sol.append(open(path_string + "/" + file, "r"))
	
	copy = open(path_string + "/copy.hpp", "w")
	tests = open(path_string + "/tests.cpp", "r")
	
	copy.write("#ifndef COPY_HPP\n")
	copy.write("#define COPY_HPP\n")
	
	for el in files:
		copy.write('#include "../solution/{}"\n'.format(el))
	
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
	
	copy.write('#include "solution.hpp"\n')
	
	parseWriteChange(student_sol, copy, tests)
	
	copy.write("#endif")
	
	copy.close()
	for item in student_sol:
		item.close()
	tests.close()



if __name__ == "__main__":
	if len(sys.argv) == 2:
		cmake(str(sys.argv[1]))
	else:
		codeexpert()
