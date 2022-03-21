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

const info_export = {
  
  	cxweb: "8.1.1",
  	cxrun: "3.12.2"
};

const permissions = {
    
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
};

const env = {
	
	slug: "generic-1",
    lastUpdateCheck: date,
    replacedBy: "generic-2"
};

async function main( ) {
	
	await fs.rmdir( export_path, { recursive: true });
	await fs.mkdir( export_path );

	await assemble_project( "solution" );
	await assemble_project( "template" );

	await fs.writeFile( `${ export_path }/export.json`, JSON.stringify( info_export, null, 2 ));

	const task_id = uid( );
	const solution_project_id = uid( );
	const template_project_id = uid( );

	const [ solution_files, solution_root_key ] = await get_dir_description( `${ export_path }/solution`, solution_project_id, true );

	//the template project files are not linked with the solution files (we do not make use of originKey)
	const [ template_files, template_root_key ] = await get_dir_description( `${ export_path }/template`, template_project_id, true );

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
			env
		},

		template: {

			_id: template_project_id,
			lastUpdate: date,
			rootDirKey: template_root_key,
			name: `${ chapter } - ${ name } - Student Template`,
			useCase: "studentTemplate",
			envVars: { },
			taskId: task_id,
			env
		},

		solutionFiles: solution_files,
		templateFiles: template_files
	};

	await fs.writeFile( `${ export_path }/task.json`, JSON.stringify( info_task, null, 2 ));
}

async function assemble_project( name ) {

	await copy_dir( `${ assignment_path }/${ name }`, `${ export_path }` );
	await copy_dir( `${ testing_path }/.`, `${ export_path }/${ name }` );
}

function copy_dir( source, target ) {

	return new Promise( resolve => {

		exec( `cp -R "${ source }" "${ target }"`, resolve );
	});
}

async function get_dir_description( path, project_id, is_top_level = false ) {

	var list = [ ];
	const children = await fs.readdir( path, { withFileTypes: true });

	const info = {

		_id: uid( ),
		key: uid( ),
		type: "inode/directory",
		version: 1 + children.length,
		current: true,
		name: is_top_level ? "." : basename( path ),
		userId: user_id,
		projectId: project_id,
		permissions,
		createdAt: date,
		children: [ ]
	};

	await Promise.all( children.map( async child => {

		const child_path = `${ path }/${ child.name }`;
		if( child.isDirectory( )) {

			const [ description, key ] = await get_dir_description( child_path );
			list = [ ...list, ...description ];
			info.children.push( key );	
		}
		else {

			const description = await get_file_description( child_path );
			list.push( description );
			info.children.push( description.key );
		}
	}));

	list.push( info );
	return [ list, info.key ];
}

async function get_file_description( path, project_id ) {

	const status = await fs.stat( path );
	const id = uid( );
	const info = {

		_id: id,
		type: await mime_type( path ),
		permissions,
		size: status.size,
		fileId: id,
		userId: user_id,
		projectId: project_id,
		version: 1,
		current: true,
		key: id,
		createdAt: date,
		name: basename( path )
	};

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