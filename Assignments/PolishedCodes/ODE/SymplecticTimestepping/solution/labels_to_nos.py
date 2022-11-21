# Python script to replace latex labels with the equation No.s
# Requires the file NCSEFL_Prbrefs.aux

import os
from os import listdir
import shutil

headerfiles = [ ff for ff in listdir("./") if ff[-3:] == "hpp"]
refs = {}

if not os.path.exists("modified_files"):
    os.mkdir("modified_files")

for hpp in headerfiles:
    shutil.copyfile(hpp,"./modified_files/"+hpp)
    with open(hpp,"r") as f:
        for line in f:
            if line.find(r'\ref{') > 0 and line.find('}') > 0:
                start = line.find(r"\ref{")
                end = line.find("}")
                refs[line[start:end+1]] = []

fname = "NCSEFL_Prbrefs.aux"

with open (fname,"r") as f:
    for line in f:
        for key in refs.keys():
            if line.find(key[5:]) > -1:
                print("Found key " + key + " in line "+ line)
                eqno = line[line.find("("):line.find(")")+1]
                refs[key] = eqno

for hpp in headerfiles:
    with open("./modified_files/"+hpp,"r") as f:
        content = ""
        for line in f:
            if line.find(r'\ref{') > 0 and line.find('}') > 0:
                start = line.find(r"\ref{")
                end = line.find("}")
                lol = line[start:end+1]
                print("Ref " + lol + " present in line \n")
                line = line.replace(lol, refs[lol])
                print("Fixed line \n " + line)
                content = content + line + "\n"
            else:
                content = content+line
    writefile = open("./modified_files/"+hpp,"w")
    writefile.write(content)


print(refs)
