# This script will format every .cpp or .hpp file,
# also in subdirectories of the directory from which it is called.

for cpp_file in $(find . -name '*.cpp' -o -name '*.hpp')
do
	# does inplace formatting using Google style guidelines
	clang-format -i --verbose --style=Google $cpp_file
done
