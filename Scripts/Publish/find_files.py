import os
import re
import sys
os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/../..")
prefix = " " + sys.argv[1] + " " if len(sys.argv)>1 else " "

#
# Configuration
#
rootdir = './'
#
# Parse acl file
#
acl = []
with open('Scripts/Publish/acl', 'r') as file:
  acl_raw = file.read().strip()
  acl_raw = acl_raw.split("\n")
  acl_raw = filter(None, acl_raw)
  # integrity check
  if len(filter(lambda acr: (acr[0]!="+" and acr[0]!="-" and acr[0] != "#"), acl_raw)) > 0:
    raise "Invalid syntax in acl file (lines must begin with + or -)"
  # remove commentary
  acl = filter(lambda acr: acr[0] != "#", acl_raw)

#
# Find all files that need to be deleted
#
files_to_be_deleted = []
# iterate over all access control rules
for acr in acl:
  # if we have a deletion rule (e.g. `-path/`) we
  #  add all files that match this rule to the
  #  files_to_be_deleted array
  if acr[0] == "-":
    path_matcher = re.compile("^" + acr[1:])
    
    # iterate over all files
    for subdir, dirs, files in os.walk(rootdir):
      for file in files:
        path = os.path.join(subdir, file)
	assert(path[0:2]=="./")
        path = path[2:]
        if path_matcher.match(path):
          files_to_be_deleted.append(path)
  elif acr[0] == "+":
    # build regex from acr
    if acr[-2:] == "/*":
      path_matcher = re.compile("^" + acr[1:-2] + "/[^/]+$")
    else:
      path_matcher = re.compile("^" + acr[1:])
    # iterate over all files
    for subdir, dirs, files in os.walk(rootdir):
      for file in files:
        path = os.path.join(subdir, file)
	assert(path[0:2]=="./")
        path = path[2:]
        if path_matcher.match(path):
          files_to_be_deleted.remove(path)
  else:
    raise "unexpected exception"

print(prefix[1:] + ("; " + prefix).join(files_to_be_deleted))
#for subdir, dirs, files in os.walk(rootdir):
#    for file in files:
#        print os.path.join(subdir, file)
