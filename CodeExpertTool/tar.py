import tarfile
import os.path

source_dir = "/Users/heinrich/repositories/NumCSE/CodeExpertTool/taskExport_Introduction_MatrixBlocks"
output_filename = "/Users/heinrich/repositories/NumCSE/CodeExpertTool/taskExport_Introduction_MatrixBlocks.tar"

with tarfile.open( output_filename, "w:gz" ) as tar:
    tar.add( source_dir, arcname = os.path.sep )