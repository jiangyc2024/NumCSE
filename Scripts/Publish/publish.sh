#!/bin/bash
#
# Variables
#
#set -o xtrace

EXIT_ON_ERROR=true
GIT_TOP_LEVEL="$(git rev-parse --show-toplevel)"
ORIGIN_BRANCH=master
CURRENT_BRANCH=$(git rev-parse --symbolic-full-name --abbrev-ref HEAD)
# runtime vars
_CLEANUP_PUBLIC_STASH=false
_CLEANUP_PUBLIC=false
#
# Utils
#
# color helpers
function heading () {
	printf '\e[48;1;4m%s\e[0m \n' "$1"
}

function notice () {
	printf '\e[0;32m%s\e[0m \n' "$1"
}

function error () {
	printf '\e[41m%s\e[0m \n' "$1"
	exit;
}

function warn () {
	printf '\e[48;5;208m%s\e[0m \n' "$1"
}


function fancy_heading () {
	#color_grad 16 21
	printf "%-60s" "$1"
	#color_grad 21 16
	printf "\n"
}

function fancy_sep () {
	local colored_ws=$(printf "%60s" " " | sed -e 's/./\\e[48;5;21m \\e[0m/g');

	#color_grad 16 21
	#printf "$colored_ws"
	#color_grad 21 16
	printf "\n"
}

# error handler
function error_handler () {
	local signal=$?
   	if $EXIT_ON_ERROR; then
		if $_CLEANUP_PUBLIC; then
		  git rebase --abort
		fi
        git checkout $CURRENT_BRANCH
		if $_CLEANUP_PUBLIC_STASH; then
		  git branch -D public_stash
		fi 
		# exit on error with error message
		[ "$signal" -eq 0 ] || error "uncaught exception occured"
	else
		[ "$signal" -eq 0 ] || warn "uncaught exception occured. continuing"
	fi
}

trap error_handler ERR # use 'trap - ERR' to remove
trap error_handler INT

#
#
#
# heading
fancy_sep;
fancy_heading "      ____________________   ___"
fancy_heading "     /  ________   ___   /__/  /"
fancy_heading "    /  _____/  /  /  /  ___   /"
fancy_heading "   /_______/  /__/  /__/  /__/"
fancy_heading "   Eidgenoessische Technische Hochschule Zuerich"
fancy_heading "   Swiss Federal Institute of Technology Zurich"
fancy_heading "------------------------------------------------------------"
fancy_heading "   Numerical Methods for CSE"
fancy_heading "    (Autumn 2016, 401-0663-00L)"
fancy_heading "------------------------------------------------------------"
fancy_heading " This script is used to publish the lecture repository"
fancy_heading " without exposing unpublished exercises or solutions."
fancy_heading ""
fancy_heading " It works as follows:"
fancy_heading "  1. Integrity checks"
fancy_heading "     Check that we are on the correct branch and have a"
fancy_heading "     clean working directory."
fancy_heading "  2. Create a branch public_stash."
fancy_heading "  3. Copy the history from the original branch into"
fancy_heading "     public_stash while removing unpublished material"
fancy_heading "  4. Rebase all commits from public_stash into public"
fancy_heading "  5. Push the changes from public to origin_public"
fancy_heading "------------------------------------------------------------"
fancy_heading " Contact"
fancy_heading "   Till Ehrengruber"
fancy_heading "   tille@student.ethz.ch"
fancy_sep;

#
# Integrity checks
#
notice "Running integrity checks"
# check that we are on the origin branch
if [ "$CURRENT_BRANCH" != "$ORIGIN_BRANCH" ]; then
  error "Fatal error. Wrong branch, please checkout $ORIGIN_BRANCH"
fi
# see http://stackoverflow.com/questions/2657935/checking-for-a-dirty-index-or-untracked-files-with-git
number_of_uncommited_files=$(git status --porcelain 2>/dev/null| egrep "^(M| M)" | wc -l)
# check the index is clean
if [ "$number_of_uncommited_files" -ne "0" ]; then
    error "Working directory is dirty. Please commit your changes or use \`git stash\` in combination with \`git stash pop\`"
fi
# check that the branch public_stash does not exist
trap - ERR # remove error handler temporary
git rev-parse --verify public_stash &> /dev/null
if [ $? -eq 0 ]; then
  error "Fatal error. Branch publish_stash already exists."
fi
trap error_handler ERR 

#
# Publish
#
notice "Remove all unpublished stuff from $ORIGIN_BRANCH into public_stash"
# go to the top level of the repo
pushd "$GIT_TOP_LEVEL" > /dev/null
# create a new branch
git branch public_stash && _CLEANUP_PUBLIC_STASH=true
# checkout new branch public stash
git checkout public_stash
# get all files to be removed and
# remove all solutions from this branch
git filter-branch -f --index-filter "$(python $GIT_TOP_LEVEL/Scripts/Publish/find_files.py 'git rm --cached --ignore-unmatch')" HEAD >> log
# check that the public branch exists
trap - ERR # remove error handler temporary
git rev-parse --verify public &> /dev/null
status=$?
trap error_handler ERR 
if [ $status -ne 0 ]; then
	echo "Creating branch public"
	git branch public
	echo "Adding remote origin_public pointing to the NumCSEStudents repo (public -> origin_public/master)"
    EXIT_ON_ERROR=false
	git remote add origin_public git@gitlab.math.ethz.ch:NumCSE/NumCSEStudent.git
	EXIT_ON_ERROR=true
fi
notice "Take all commits from public_stash and rebase them onto public"
# checkout the public branch
git checkout public
# incorporate all changes
_CLEANUP_PUBLIC=true
git rebase public_stash
_CLEANUP_PUBLIC=false
# push
# first check if public and public_stash share the same history
if [ $(git rev-parse public) != $(git rev-parse public_stash) ]; then
  warn "Your branches have diverged!"
  warn "Did you accidently exclude something from the repo in \`acl\` that was already published?"
  warn "Do you still want to push? (use at your own risk)"
  read -r -p "[y/N]: " response
	case $response in
		[yY][eE][sS]|[yY])
            echo "Pushing changes"
			git push -f -u origin_public public:master
			;;
	esac
else
  echo "Pushing changes..."
  git push -f -u origin_public public:master
fi
# remove public_stash
git branch -D public_stash
_CLEANUP_PUBLIC_STASH=false
# checkout master
git checkout master
notice "Finished."

popd > /dev/null
