import sys

wdl = sys.argv[1]
fd = open (wdl)
for line in fd:
	if "docker:" not in line:
		print line,