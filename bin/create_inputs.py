# Python Script to Prepare WDL Inputs
import sys

def read_config():
	config_map = {}
	fd = open("conf/config.txt")
	for line in fd:
		line = line.strip()
		parts = line.split("=")
		config_map[parts[0].strip()] = parts[1].strip()
	return config_map

json = sys.argv[1]
environment = sys.argv[2]
config_map = read_config()
google_bucket = config_map["gs_bucket"]

env_base_url = ""
if environment == "google":
	env_base_url = google_bucket
fd = open (json)
for line in fd:
	line = line.replace("${base_dir}/", env_base_url)
	print line,
