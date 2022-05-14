/**

Authors: 

	Heinrich Grattenthaler (heinrich.grattenthaler@sam.math.ethz.ch)
	Pascal Müller (pamuelle@ethz.ch)

Disambiguation:
	
	cx = CodeExpert (expert.ethz.ch)
	task = cx terminology for assignment
	project = either solution or template, cx terminology for the two subfolders of an assignment
	environment = the Docker container + some scripts + configuration used by cx to run a project
	key = unique identifiers used for files/directories to indicate the structure of a directory

Usage: 

	node export.js <repo_root> <assignment_path_rel> [<default_file>] [--keep]

This script reverse-engineers the cx task export format for automating part of the process of 
uploading assignments to the cx platform. It produces a .tar archive ready for import in cx. 
the cx template project files are typically linked to their corresponding solution files (by using
the keys of the solution files) so we make an effort to do the same. Otherwise the template project
will be marked faulty after import and cannot be published. See details below.

**/

const fs = require( "fs/promises" );
const { exec } = require( "child_process" );
const { basename } = require( "path" );

//absolute path to root of the repository
const repo_root = process.argv[ 2 ];

//path to the assignment directory, relative from repository root
const assignment_path_rel = process.argv[ 3 ]

//the basename of the default student work file (optional, only needed if the auto-detection fails 
//and you get an error)
const default_file = process.argv[ 4 ];

//whether to keep the export directory (unzipped)
const keep = process.argv.includes( "--keep" );

//display name of the task
const display_name = basename( assignment_path_rel );

//path/name used for the export directory/archive, relative to working directory
const export_path = "./taskExport";

//creation date used for all file meta data
const date = ( new Date( )).toISOString( );

//some random cx user id, likely not important
const user_id = "au4g3HpKS9AbksWyF";

//for name-based creation of links between solution and template
const key_cache = { };

//for finding the default file
const default_cache = { };

//version info supplied to cx about the imported data format
const info_export = {
  
  	cxweb: "8.1.1",
  	cxrun: "3.12.2"
};

//the strict, default permissions for all files
const permissions_strict = _ => ({
    
  phase: [
    "interactive",
    "submission"
  ],
  read: [
    "admin",
    "assistant"
  ],
  write: [
    "admin"
  ]
});

//the cx environment info used for this task
const env = {
	
	slug: "generic-1",
	lastUpdateCheck: date,
  replacedBy: "generic-2"
};

//detection logic for whether this is the default student work file
async function is_default_file( path ) {

	return ( default_file && basename( path ) == basename( default_file )) ? true : await contains_todo( path );
}

//detection logic for whether a file is marked as writable by students in cx
async function is_editable( path ) {

	if( await is_default_file( path )) return true;
	if( await contains_todo( path )) return true;
	if( path.includes( "/written_solution.md" ) || path.includes( "main.cpp" )) return true;
}

//helper to determine whether this is a c++ file containing work for students
async function contains_todo( path ) {

	if( path.includes( ".hpp" ) && ! path.endsWith( "solution.hpp" )) {

		const code = await fs.readFile( path );
		if( code.includes( "TODO" ) || code.includes( "TO DO" )) {

			return true;
		}
	}
}

//detection logic for whether a file is marked as readable by students in cx
async function is_readable( path ) {

	if( await is_editable( path )) return true;
	if(( path.endsWith( ".h" ) || path.endsWith( ".hpp" )) && ! path.endsWith( "solution.hpp" )) return true;
}

//paths of c++ files that are added to the projects by default, relative to the repository root
const default_includes = [ 

	"MatplotlibC++/matplotlibcpp.h",
	"Utils/timer.h"
];

//-------------------- BUSINESS LOGIC --------------------//

const solution_project_id = uid( );
const template_project_id = uid( );

async function main( ) {
	
	//clean directory if existing
	await fs.rm( export_path, { recursive: true, force: true });
	await fs.mkdir( export_path );

	//copy necessary files into both projects
	await assemble_project( "solution" );
	await assemble_project( "template" );

	await fs.writeFile( `${ export_path }/export.json`, JSON.stringify( info_export, null, 2 ));

	//obtain the file structure in cx serialization format and keys of the project directories for both projects
	//it is vital that the solution project is processed first, since the linking from template to solution files relies on it
	const [ solution_files, solution_root_dir_key ] = await get_dir_description( `${ export_path }/solution`, solution_project_id, "." );
	const [ template_files, template_root_dir_key ] = await get_dir_description( `${ export_path }/template`, template_project_id, "." );
	
	//build the main info file for all cx task data
	const task_id = uid( );
	const info_task = { 

		_id: task_id,
		name: display_name,
		kind: "code",
		masterSolutionId: solution_project_id,
		studentTemplateId: template_project_id,
		
		solution: {

			_id: solution_project_id,
			lastUpdate: date,
			rootDirKey: solution_root_dir_key,
			name: `${ display_name } - Master Solution`,
			useCase: "masterSolution",
			envVars: { },
			taskId: task_id,
			env,
			defaultSelectionKey: default_cache[ solution_project_id ]
		},

		template: {

			_id: template_project_id,
			lastUpdate: date,
			rootDirKey: template_root_dir_key,
			name: `${ display_name } - Student Template`,
			useCase: "studentTemplate",
			envVars: { },
			taskId: task_id,
			env,
			defaultSelectionKey: default_cache[ template_project_id ]
		},

		solutionFiles: solution_files,
		templateFiles: template_files
	};

	await fs.writeFile( `${ export_path }/task.json`, JSON.stringify( info_task, null, 2 ));
	await zip( export_path );

	//comment out for debugging the produced archive
	if( ! keep ) await fs.rm( export_path, { recursive: true });
}

async function assemble_project( name ) {

	//copy the project directory into our task
	await copy_dir( `${ repo_root }/${ assignment_path_rel }/${ name }`, `${ export_path }` );
	
	const project_path = `${ export_path }/${ name }`;
	const testing_path = `${ repo_root }/Testing`;

	//files necessary to conduct the tests in cx
	await copy_dir( `${ testing_path }/scripts`, project_path );
	await copy_file( `${ testing_path }/conf.yml`, project_path );
	await copy_file( `${ testing_path }/doctest.h`, project_path );
	await copy_file( `${ testing_path }/copy_and_tweak.py`, project_path );

	//find the master solution by finding the default student work file in the solution project. 
	//solution.hpp is needed for the tests
	const solution_path = await find_default_file( `${ export_path }/solution`, true );
	await copy_file( solution_path, `${ project_path }/solution.hpp` );

	//just a parallel for-loop, copy default library files
	await Promise.all( default_includes.map( async file => {

		await copy_file( `${ repo_root }/${ file }`, project_path );
	}));
}

//async copy directory wrapper
function copy_dir( source, target ) {

	return new Promise( resolve => {

		exec( `cp -R "${ source }" "${ target }"`, resolve );
	});
}

//async copy file wrapper
function copy_file( source, target ) {

	return new Promise( resolve => {

		exec( `cp "${ source }" "${ target }"`, resolve );
	});
}

//required because cx only accepts .tar archives with a specific internal naming convention
function zip( dir ) {

	return new Promise( resolve => {

		exec( `python tar.py "${ dir }"`, resolve );
	});
}

//tries to find the default student work file in a directory
async function find_default_file( path, is_top_level = false ) {

	var result = undefined;
	const children = await fs.readdir( path, { withFileTypes: true });
	
	//clean up and exit with error
	const fail = async _ => {

		await fs.rmdir( export_path, { recursive: true });
		console.error( "Error: Auto-detection of default file failed. Please specify it explicitly in the arguments" );
		process.exit( 1 );
	};

	//just a parallel for-loop
	await Promise.all( children.map( async child => {

		const child_path = `${ path }/${ child.name }`;

		if( child.isDirectory( )) {

			const sub_result = await find_default_file( child_path );
			if( result && sub_result ) await fail( );
			if( sub_result ) result = sub_result;
		}
		else {

			if( await is_default_file( child_path )) {

				//duplicate found
				if( result ) await fail( );
			 	result = child_path;
			}
		}
	}));

	//none found in the whole recursive search
	if( ! result && is_top_level ) await fail( );
	return result;
}

//builds a directory description in cx format, cx_path is the project-relative path used by the
//internal cx format
async function get_dir_description( path, project_id, cx_path ) {

	//accumulates file system node information in this directory (recursively) in a flat list
	const list = [ ];
	const children = await fs.readdir( path, { withFileTypes: true });

	//info about this directory
	const info = {

		_id: uid( ),
		key: uid( ),
		type: "inode/directory",
		version: 1 + children.length,
		current: true,
		name: basename( cx_path ),
		userId: user_id,
		projectId: project_id,
		permissions: permissions_strict( ),
		createdAt: date,
		children: [ ]
	};

	//used for linking the template project directories to the solution project directories by cx_path
	const cached = key_cache[ cx_path ];
	
	if( cached ) {

		//there exists a directory with the same name in the solution project, so create a link
		info.originKey = cached.key;
		info.originVersion = cached.version;
	}
	else {

		//this is the first directory with this cx_path so create an entry for the template to link to
		key_cache[ cx_path ] = info;
	}

	//just a parallel for-loop
	await Promise.all( children.map( async child => {

		//file system path and internal cx path
		const child_path = `${ path }/${ child.name }`;
		const child_cx_path = `${ cx_path }/${ child.name }`;

		if( child.isDirectory( )) {

			//get node information from child directory and its key
			const [ description, key ] = await get_dir_description( child_path, project_id, child_cx_path );
			info.children.push( key );

			//append all grand-children's descriptions to the global list
			list.push( ...description );
		}
		else {

			//get node information from child file
			const description = await get_file_description( child_path, project_id, child_cx_path );
			info.children.push( description.key );

			//append child description to the global list
			list.push( description );
		}
	}));

	//push own node info to the global list
	list.push( info );
	return [ list, info.key ];
}

async function get_file_description( path, project_id, cx_path ) {

	const id = uid( );
	
	//info about this file
	const info = {

		_id: id,
		type: await mime_type( path ),
		permissions: permissions_strict( ),
		size: ( await fs.stat( path )).size,
		fileId: id,
		userId: user_id,
		projectId: project_id,
		version: 1,
		current: true,
		key: id,
		createdAt: date,
		name: basename( cx_path )
	};

	//used for linking the template project files to the solution project files by cx_path
	const cached = key_cache[ cx_path ];

	if( cached ) {

		//there exists a file with the same name in the solution project, so create a link
		info.originKey = cached.key;
		info.originVersion = cached.version;

		//additionally indicate that you wish to ignore the file content in the template and use the 
		//solution file directly
		info.fileId = cached.key;
	}
	else {

		//this is the first file with this cx_path so create an entry for the template to link to
		key_cache[ cx_path ] = info;
	}

	//save the info about the default file for each project, so the projects can point to this file 
	//as default file
	if( await is_default_file( path )) default_cache[ project_id ] = info.key;
	
	if( await is_editable( path )) {

		//add permissions to make this file writable (for students) in cx
		info.permissions.write.push( "assistant", "student", "everyone" );
		
		//if this is the template, do not use the solution file, but keep a divergent version from the 
		//solution. the link is still functional. cx will tell you that about said diversion when you 
		//click the file in the GUI, but this is intended behavior
		if( project_id == template_project_id ) info.fileId = info.key;
	}

	//add permissions to make this file readable (for students) in cx
	if( await is_readable( path )) info.permissions.read.push( "student", "everyone" );

	return info;
}

//async wrapper for detecting the mime type of a file. not sure whether this matters beyond a
//cosmetic level. the output is not completely aligned with the mime types that cx detects when
//you export a task
async function mime_type( path ) {

	return new Promise( resolve => {

		exec( `file --mime-type -b "${ path }"`, ( _, stdout ) => resolve( stdout.replaceAll( "\n", "" )));
	});
}

//generate a unique identifier with the same format as cx's uids
function uid( ) {

	var str = "";
	const chars = "ABCDEFGHIHJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

	for( var i = 0; i < 17; ++ i ) {

		const index = Math.floor( Math.random( ) * chars.length );
		str += chars[ index ];
	}

	return str;
}

main( );