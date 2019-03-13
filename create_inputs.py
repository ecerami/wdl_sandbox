import sys

json = sys.argv[1]
environment = sys.argv[2]
env_base_url = ""
if environment == "google":
	env_base_url = "gs://cromwellbucket4221/"
fd = open (json)
for line in fd:
	line = line.replace("${base_dir}/", env_base_url)
	print line,