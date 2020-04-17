import glob
import re
import sys

#
# CURRENT ISSUES:
#   - Does not spot student's helper functions.
#   - The cmake version can only be called from the folder's local cmake file and its build folder.
#   - Can only handle one file that is written by the student.
#   - probably more issues
#

def parseWriteChange(student_sol, copy):
	"""
	Checks the c++ file student_sol for function definitions and changes the name of those to test_*.
	"""
	prefix = "test_"
	
	for line in student_sol:
		# parse file for function names and change them
		parse1 = re.match("([a-zA-Z0-9_\<\>\:]+)(\s+)([a-zA-Z0-9_]*)(\()(.*)", line)
		parse2 = re.match("([a-zA-Z0-9_]*)(\()(.*)", line)
		if parse1 and not re.match("using|template|namespace", line):
			copy.write(parse1.group(1) + parse1.group(2) + prefix + parse1.group(3)\
			+ parse1.group(4) + parse1.group(5))
		elif parse2 and not re.match("using|template|namespace", line):
			copy.write(prefix + parse2.group(1) + parse2.group(2) + parse2.group(3))
		else:
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
	expected_hpps_wout_student_sol = set(["../polyfit.hpp", "../copy.hpp"])
	all_hpps = set(glob.glob("../*.hpp"))
	
	file = compareHeadersWithExpected(expected_hpps_wout_student_sol, all_hpps)
	
	student_sol = open(file, "r")
	copy = open("../copy.hpp", "w")
	
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
	
	parseWriteChange(student_sol, copy)
	
	if no_include_guards:
		copy.write("#endif")
	
	copy.close()
	student_sol.close()

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
	
	parseWriteChange(student_sol, copy)
	
	if no_include_guards:
		copy.write("#endif")
	
	copy.close()
	student_sol.close()



if __name__ == "__main__":
	if len(sys.argv) == 2:
		if str(sys.argv[1]) == "cmake":
			cmake()
		else:
			raise SystemExit("Wrong usage: either 'cmake' as argument to run cmake version or \
			nothing implying the CodeExpert version.")
	else:
		codeexpert()
