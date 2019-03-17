# Python Script to Strip Docker attributes from a WDL File
import sys
wdl = sys.argv[1]
fd = open (wdl)
for line in fd:
	if "docker:" not in line:
		print line,