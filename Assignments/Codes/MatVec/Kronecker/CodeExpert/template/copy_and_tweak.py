import glob
import re
import sys

if __name__ == "__main__":
	
	arg = ""
	cmake_helper = ""
	if len(sys.argv) == 2:
		arg = str(sys.argv[1])
		if arg == "cmake":
			cmake_helper = "../"
	
	expected_hpps_wout_student_sol = set([cmake_helper + "polyfit.hpp", cmake_helper + "copy.hpp"])
	
	all_hpps = set(glob.glob(cmake_helper + "*.hpp"))
	diff = all_hpps.difference(expected_hpps_wout_student_sol)
	
	# get the student header
	if len(diff) == 1:
		file = diff.pop()
	else:
		raise SystemExit("Problem finding the student's solution.")
	
	student_sol = open(file, "r")
	
	copy = open(cmake_helper + "copy.hpp", "w")
	
	# expect include guards, change them
	if(student_sol.readline().startswith("#ifndef")):
		next(student_sol)
		next(student_sol)
		no_include_guards = False
	else:
		no_include_guards = True
	copy.write("#ifndef COPY_HPP\n")
	copy.write("#define COPY_HPP\n")
	
	# include the solution
	if arg == "cmake":
	    copy.write('#include "../solution/{}"\n'.format(file[3:]))
	else:
	    copy.write('#include "../solution/{}"\n'.format(file))
	
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
	
	if no_include_guards:
		copy.write("#endif")

	copy.close()
	student_sol.close()
