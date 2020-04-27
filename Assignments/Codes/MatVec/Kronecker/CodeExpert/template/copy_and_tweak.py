import re
import sys

#
# CURRENT ISSUES:
#   - The cmake version can only be called from the folder's local cmake file and its build folder.
#   - does not recognize inheriting classes
#   - probably more issues
#

def parseWriteChange(student_sol, copy, tests):
	"""
	Checks the c++ file student_sol for function definitions and changes the name of those to test_*.
	"""
	suffix = "_TEST"
	
	# get all function signatures to be tested from the test file
	function_signatures = []
	for line in tests:
		parse = re.match('\s*TEST_CASE\("([\S\s]*)(\s*)(\S*)"\s\*\sdoctest::description.*', line)
		if parse:
			function_signatures.append(parse)
	
	for file in student_sol:
		# expect include guards, change them
		if file.readline().startswith("#ifndef"):
			next(file)
			next(file)
		elif file.readline().startswith("#pragma once"):
			next(file)
		for line in file:
			write = True
			for el in function_signatures:
				parse = re.match("(\s*)" + el.group(1) + el.group(2) + el.group(3) + "(\s*\(.*)", line)
				if parse:
					copy.write(parse.group(1) + el.group(1) + el.group(2) + suffix + el.group(3) + parse.group(2) + "\n")
					write = False
			parse_class = re.match("\s*(class|struct)(\s*)(\S*)(\s*)({|\s*)", line)
			if parse_class:
				copy.write(parse_class.group(1) + parse_class.group(2) + parse_class.group(3) + suffix + parse_class.group(4) + parse_class.group(5) + "\n")
				write = False
			if write and not line.startswith("#endif"):
				copy.write(line)
		copy.write("\n")
				
def searchStudentSolution(main):
	"""
	Searches the student solution using main.cpp and the exact order of inclusion in it.
	Returns a list of filenames.
	"""
	return_list = []
	other_headers = ["polyfit", "solution", "solution2", "solution3"] # other headers ending with .hpp
	for line in main:
		parse = re.match('\s*\#include\s*"(.*).hpp"', line)
		if parse and parse.group(1) not in other_headers:
			return_list.append(parse.group(1) + ".hpp")
	return return_list

def cmake():
	"""
	Will perform the copy and tweak operation if called from a cmake build system inside a build folder.
	The tweaked copy is named 'copy.hpp'.
	Will also put header guards into the file.
	"""
	main = open("../main.cpp", "r")
	files = searchStudentSolution(main)
	main.close()
	student_sol = []
	for file in files:
		student_sol.append(open("../" + file, "r"))
	
	copy = open("../copy.hpp", "w")
	tests = open("../tests.cpp", "r")
	
	copy.write("#ifndef COPY_HPP\n")
	copy.write("#define COPY_HPP\n")
	
	copy.write('#include "../solution/{}"\n'.format(file[3:]))
	
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
		if str(sys.argv[1]) == "cmake":
			cmake()
		else:
			raise SystemExit("Wrong usage: either 'cmake' as argument to run cmake version or \
			nothing implying the CodeExpert version.")
	else:
		codeexpert()
