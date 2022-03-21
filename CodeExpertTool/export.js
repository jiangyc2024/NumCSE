const fs = require( "fs/promises" );
const { exec } = require( "child_process" );
const { basename } = require( "path" );
const assignment_path = process.argv[ 2 ]
const testing_path = process.argv[ 3 ];	
const build_path = process.argv[ 4 ];
const chapter = "Introduction"
const name = "MatrixBlocks"
const date = ( new Date( )).toISOString( );
const export_path = `${ build_path }/${[ "taskExport", chapter, name ].join( "_" )}`;
const user_id = "au4g3HpKS9AbksWyF";
const solution_project_id = uid( );
const template_project_id = uid( );

//for name-based memorizing to create links between solution and template
const key_cache = { };

//for finding the default file
const default_cache = { };

const info_export = {
  
  	cxweb: "8.1.1",
  	cxrun: "3.12.2"
};

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

const env = {
	
	slug: "generic-1",
    lastUpdateCheck: date,
    replacedBy: "generic-2"
};

async function is_default_file( path ) {

	if( path.endsWith( ".hpp" ) && ! path.endsWith( "solution.hpp" )) {

		const code = await fs.readFile( path );
		if( code.includes( "TODO:" ) && code.includes( "SAM_LISTING" )) {

			return true;
		}
	}
}

async function is_editable( path ) {

	if( await is_default_file( path )) return true;
	if( path.includes( "/written_solution.md" ) || path.includes( "main.cpp" )) return true;
}

async function is_readable( path ) {

	if( await is_editable( path )) return true;
	if(( path.endsWith( ".h" ) || path.endsWith( ".hpp" )) && ! path.endsWith( "solution.hpp" )) return true;
}

//-------------------- BUSINESS LOGIC --------------------//

async function main( ) {
	
	await fs.rmdir( export_path, { recursive: true });
	await fs.mkdir( export_path );

	await assemble_project( "solution" );
	await assemble_project( "template" );

	await fs.writeFile( `${ export_path }/export.json`, JSON.stringify( info_export, null, 2 ));

	const [ solution_files, solution_root_key ] = await get_dir_description( `${ export_path }/solution`, solution_project_id, "." );
	const [ template_files, template_root_key ] = await get_dir_description( `${ export_path }/template`, template_project_id, "." );
	
	const task_id = uid( );
	const info_task = { 

		_id: task_id,
		name: `${ chapter } - ${ name }`,
		kind: "code",
		masterSolutionId: solution_project_id,
		studentTemplateId: template_project_id,
		
		solution: {

			_id: solution_project_id,
			lastUpdate: date,
			rootDirKey: solution_root_key,
			name: `${ chapter } - ${ name } - Master Solution`,
			useCase: "masterSolution",
			envVars: { },
			taskId: task_id,
			env,
			defaultSelectionKey: default_cache[ solution_project_id ]
		},

		template: {

			_id: template_project_id,
			lastUpdate: date,
			rootDirKey: template_root_key,
			name: `${ chapter } - ${ name } - Student Template`,
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
}

async function assemble_project( name ) {

	await copy_dir( `${ assignment_path }/${ name }`, `${ export_path }` );
	
	const project_path = `${ export_path }/${ name }`;
	await copy_dir( `${ testing_path }/scripts`, project_path );
	await copy_file( `${ testing_path }/conf.yml`, project_path );
	await copy_file( `${ testing_path }/doctest.h`, project_path );
	await copy_file( `${ testing_path }/copy_and_tweak.py`, project_path );

	const solution_path = await find_default_file( `${ export_path }/solution` );
	await copy_file( solution_path, `${ project_path }/solution.hpp` );
}

function copy_dir( source, target ) {

	return new Promise( resolve => {

		exec( `cp -R "${ source }" "${ target }"`, resolve );
	});
}

function copy_file( source, target ) {

	return new Promise( resolve => {

		exec( `cp "${ source }" "${ target }"`, resolve );
	});
}

async function find_default_file( path ) {

	var result = undefined;
	const children = await fs.readdir( path, { withFileTypes: true });

	await Promise.all( children.map( async child => {

		const child_path = `${ path }/${ child.name }`;

		if( child.isDirectory( )) {

			result ??= await find_default_file( child_path );
		}
		else {

			if( await is_default_file( child_path )) result = child_path;
		}
	}));

	return result;
}

async function get_dir_description( path, project_id, dir_name ) {

	var list = [ ];
	const children = await fs.readdir( path, { withFileTypes: true });

	const info = {

		_id: uid( ),
		key: uid( ),
		type: "inode/directory",
		version: 1 + children.length,
		current: true,
		name: basename( dir_name ),
		userId: user_id,
		projectId: project_id,
		permissions: permissions_strict( ),
		createdAt: date,
		children: [ ]
	};

	const cached = key_cache[ dir_name ];
	
	if( cached ) {

		info.originKey = cached.key;
		info.originVersion = cached.version;
	}
	else {

		key_cache[ dir_name ] = info;
	}

	await Promise.all( children.map( async child => {

		const child_path = `${ path }/${ child.name }`;
		const child_name = `${ dir_name }/${ child.name }`;

		if( child.isDirectory( )) {

			const [ description, key ] = await get_dir_description( child_path, project_id, child_name );
			list = [ ...list, ...description ];
			info.children.push( key );	
		}
		else {

			const description = await get_file_description( child_path, project_id, child_name );
			list.push( description );
			info.children.push( description.key );
		}
	}));

	list.push( info );
	return [ list, info.key ];
}

async function get_file_description( path, project_id, file_name ) {

	const status = await fs.stat( path );
	const id = uid( );
	
	const info = {

		_id: id,
		type: await mime_type( path ),
		permissions: permissions_strict( ),
		size: status.size,
		fileId: id,
		userId: user_id,
		projectId: project_id,
		version: 1,
		current: true,
		key: id,
		createdAt: date,
		name: basename( file_name )
	};

	const cached = key_cache[ file_name ];

	if( cached ) {

		info.originKey = cached.key;
		info.originVersion = cached.version;
		info.fileId = cached.key;
	}
	else {

		key_cache[ file_name ] = info;
	}

	if( await is_default_file( path )) default_cache[ project_id ] = info.key;
	if( await is_editable( path )) {

		info.permissions.write.push( "assistant", "student", "everyone" );
		
		//will diverge from solution/origin, but keep the connection
		if( project_id == template_project_id ) info.fileId = info.key;
	}

	if( await is_readable( path )) info.permissions.read.push( "student", "everyone" );

	return info;
}

async function mime_type( path ) {

	return new Promise( resolve => {

		exec( `file --mime-type -b "${ path }"`, ( _, stdout ) => resolve( stdout.replaceAll( "\n", "" )));
	});
}

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