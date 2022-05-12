import tarfile
import os.path
import sys

source_dir = sys.argv[ 1 ]
output_filename = os.path.basename( source_dir ) + ".tar"

with tarfile.open( output_filename, "w:gz" ) as tar:
    
    tar.add( source_dir, arcname = os.path.sep )