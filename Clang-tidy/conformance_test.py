"""
Runs clang-tidy on NumCSE codes and filters out warnings from selected third-party includes. 
This script requires exported compile_commands.json from the cmake generation step.
Usage: python conformance_test.py [-f <file> |-d <directory>] [-p <build_directory>] [--no_filter] [--format]
	-f: process only a single file (<file>)
	-d: process a whole directory (<directory>)
	neither -f nor -d: process the working directory
	-p: specify the build directory (<build_directory>). This is the working directory by default.
	--no_filter: do not filter out uninteresting warnings
	--format: additionally run clang-format on all processed files
"""

import os, sys, io, shutil, math, re, subprocess

build_path = os.path.abspath( "." )
delim = '\d+:\d+: (warning|error)'
no_filter = False
do_format = False

( w,_ ) = shutil.get_terminal_size(( 50, 20 ))

def headline( title ):
	"""Print a headline"""
	print( "="*int((w-len(title))/2-1) + " " + title + " " + "="*int(math.ceil((w-len(title))/2-1 )))

def capture_cmd( cmd ):
	"""Run a subprocess and return the stdout output"""
	p = subprocess.run( cmd, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = True, check = False )
	return p.stdout.decode( 'utf-8' )

path_pat = '(?P<path>.*)(?P<delim>'+ delim +')'

def is_relevant( warning ):
	"""Check if a warning should be shown to user. If no_filter is True, every warning is shown. Errors are always displayed."""
	matches = re.search( path_pat, warning )
	if not matches:
		return False
	if no_filter:
		return True
	path = matches.group( 'path' )
	delim_match = matches.group( 'delim' )
	if 'error' in delim_match:
		return True
	return not 'MatplotlibC++' in path and not 'Eigen_install' in path	

bold_magenta = '\033[35;1m'
bold_red = '\033[31;1m'
bold = '\033[1m'
green = '\033[32m'
normal = '\033[0m'

def pretty( w ):
	"""Return a warning with some color for better readability"""
	w = w.replace( 'warning:', bold_magenta + 'warning:' + normal )
	w = w.replace( 'error:', bold_red + 'error:' + normal )
	w = w.replace( 'note:', bold + 'note:' + normal )
	w = w.replace( '^', green + '^' + normal )
	return w

def warn_filter( m ):
	"""Return a filtered version of the clang-tidy output"""
	if m == 'error': return ''
	warnings = []
	( w, m ) = parse_warning( m )
	c = 0
	while w:
		if is_relevant( w ):
			warnings.append( w )
		( w, m ) = parse_warning( m )
		c += 1
	log = "Suppressed " + str( c - len( warnings )) + " out of " + str( c )
	if( not len( warnings )): log += ' ✅'
	for w in warnings:
		log += "\n\n" + pretty( w )
	return log

pat = '(?P<warn>[^\n]*?(' + delim + ')(.|\n)*?)(?=([^\n]*?(' + delim + '))|$)'

def parse_warning( m ):
	"""Consume/parse one warning and return result and unconsumed message in typical parser fashion"""
	matches = re.search( pat, m )
	if not matches:
		return ('',m)
	warning = matches.group( 'warn' )
	left = m.replace( warning, '', 1)
	return (warning, left)

def go( item ):
	"""Process a file and print results"""
	( w,_ ) = shutil.get_terminal_size(( 50, 20 ))
	print( "-"*w )
	print( "File: {}".format( item ))
	print( )
	if do_format:
		print( "* clang-format:" )
		os.system( "clang-format -style='Google' -i '{}'".format( item ))
		print( "Done\n")
	print( "* clang-tidy:" )
	s = capture_cmd( "clang-tidy '{0}' -p {1}".format( item, build_path ))
	print( warn_filter( s ))
	print( "Done" )
	print( "-"*w )

def go_dir( path ):
	"""Process all files in a directoy"""
	os.chdir( path )
	for item in os.listdir( path ):
		fullpath = path + "/" + item
		if os.path.isfile( fullpath ) and (item.endswith( ".cpp" ) or item.endswith( ".hpp" )):
			go( item )

def main( ):
	global no_filter, do_format, build_path
	headline( "Conformance Test" )
	if "--no_filter" in sys.argv: 
		no_filter = True

	if "--format" in sys.argv:
		do_format = True

	if "-p" in sys.argv:
		i = sys.argv.index( "-p" )
		if( len( sys.argv ) >= i + 2 ):
			build_path = os.path.abspath( sys.argv[ i + 1 ])

	if len( sys.argv ) >= 3:
		if sys.argv[ 1 ] == "-f":
			go( sys.argv[ 2 ] )	
		elif sys.argv[ 1 ] == "-d": 
			go_dir( sys.argv[ 2 ] )
		else: 
			go_dir( "." )
	else:
		go_dir( "." )

	headline( "Done" )

if __name__ == '__main__':
	main( )