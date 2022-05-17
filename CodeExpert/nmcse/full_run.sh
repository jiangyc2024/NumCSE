#! /bin/bash
# 
# Usage: ./full_run.sh PolynomialInterpolation/EvaluatingDerivatives test
#        ./full_run.sh <shortened, relative assignment path> [<action>]
# 
# This is a utility wrapper around export.sh and run.sh. 
# It exports a code from our repository and then runs it in the code expert simulator, all in one stop


parent_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P)"
tool_path="$parent_dir/../tool"
$tool_path/export.sh $1 --keep;
rm -r $parent_dir/projectfiles/*;
cp -r $tool_path/taskExport/solution/* $parent_dir/projectfiles/;
$parent_dir/run.sh ${@:2};
