import re
import sys
import os
import shutil


# Read in the pyproject.toml file to get the current version
with open("pyproject.toml") as file:
    for line in file:
        if "version = " in line:
            m = re.search(R'"([A-Za-z0-9_\./\\-]*)"', line)
            ver = list(map(int, m.group(0).replace('"', '').split(".")))


# Increment according to the cmdline argument
try:
    v = sys.argv[1]
except:
    v = "micro"

if v == "micro":
    ver[2] += 1
elif v == "minor":
    ver[1] += 1
elif v == "major":
    ver[0] += 1
else:
    raise ValueError("Version to bump must be one of 'micro', 'minor', 'major'")


# Write the new version to the pyproject.toml file
shutil.copy("pyproject.toml", "pyproject.toml.tmp")
with open("pyproject.toml", mode="w") as act:
    with open("pyproject.toml.tmp") as tmp:
        for line in tmp:
            if "version = " in line:
                act.write(f'version = "{".".join(map(str, ver))}"\n')
            else:
                act.write(line)
os.remove("pyproject.toml.tmp")

