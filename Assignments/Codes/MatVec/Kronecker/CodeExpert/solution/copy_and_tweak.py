import glob
import re
import sys

#
# CURRENT ISSUES:
#   - The cmake version can only be called from the folder's local cmake file and its build folder.
#   - Can only handle one file that is written by the student.
#   - probably more issues
#

def parseWriteChange(student_sol, copy, tests):
	"""
	Checks the c++ file student_sol for function definitions and changes the name of those to test_*.
	"""
	prefix = "test_"
	
	# get all function signatures to be tested from the test file
	function_signatures = []
	for line in tests:
		parse = re.match('[\s\t]TEST_CASE\("(\S*)([\s]*)(\S*)"\s\*\sdoctest::description.*', line)
		if parse:
			function_signatures.append(parse)
	
	for line in student_sol:
		write = True
		for el in function_signatures:
			parse = re.match("([\s]*)" + el.group(1) + el.group(2) + el.group(3) + "([\s]*\(.*)", line)
			if parse:
				copy.write(parse.group(1) + el.group(1) + el.group(2) + prefix + el.group(3) + parse.group(2) + "\n")
				write = False
		if write:
			copy.write(line)
			
def compareHeadersWithExpected(expected, all):
	"""
	Compares the set expected with all and returns the file student_sol.
	"""
	diff = all.difference(expected)
	# get the student header
	if len(diff) == 1:
		return diff.pop()
	else:
		raise SystemExit("Problem finding the student's solution.")

def cmake():
	"""
	Will perform the copy and tweak operation if called from a cmake build system inside a build folder.
	The tweaked copy is named 'copy.hpp'.
	Will also put header guards into the file.
	"""
	expected_hpps_wout_student_sol = set(["../polyfit.hpp", "../copy.hpp", "../solution.hpp"])
	all_hpps = set(glob.glob("../*.hpp"))
	
	file = compareHeadersWithExpected(expected_hpps_wout_student_sol, all_hpps)
	
	student_sol = open(file, "r")
	copy = open("../copy.hpp", "w")
	tests = open("../tests.cpp", "r")
	
	# expect include guards, change them
	if(student_sol.readline().startswith("#ifndef")):
		next(student_sol)
		next(student_sol)
		no_include_guards = False
	else:
		no_include_guards = True
	copy.write("#ifndef COPY_HPP\n")
	copy.write("#define COPY_HPP\n")
	
	copy.write('#include "../solution/{}"\n'.format(file[3:]))
	
	parseWriteChange(student_sol, copy, tests)
	
	if no_include_guards:
		copy.write("#endif")
	
	copy.close()
	student_sol.close()
	tests.close()

def codeexpert():
	"""
	Will perform the copy and tweak operation if called from a script that is used by codeexpert.
	The tweaked copy is named 'copy.hpp'.
	Will also put header guards into the file.
	Puts copy.hpp in a writable folder on the server – MUST BE REMOVED AFTER COMPILING!
	"""
	expected_hpps_wout_student_sol = set(["polyfit.hpp", "solution.hpp"])
	all_hpps = set(glob.glob("*.hpp"))
	
	file = compareHeadersWithExpected(expected_hpps_wout_student_sol, all_hpps)
	
	student_sol = open(file, "r")
	copy = open("copy.hpp", "w")
	tests = open("tests.cpp", "r")
	
	# expect include guards, change them
	if(student_sol.readline().startswith("#ifndef")):
		next(student_sol)
		next(student_sol)
		no_include_guards = False
	else:
		no_include_guards = True
	copy.write("#ifndef COPY_HPP\n")
	copy.write("#define COPY_HPP\n")
	
	copy.write('#include "solution.hpp"\n')
	
	parseWriteChange(student_sol, copy, tests)
	
	if no_include_guards:
		copy.write("#endif")
	
	copy.close()
	student_sol.close()
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
