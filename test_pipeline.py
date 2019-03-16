import subprocess
import shlex
import json
import time

cromwell_base = "http://docker.for.mac.localhost:8000/api/workflows/v1"

# Run a subprocess and capture the JSON Response
def run_subprocess(command_line):
	args = shlex.split(command_line)
	sub = subprocess.check_output(args)
	parsed_json = json.loads(sub)
	return (parsed_json)

# Submit new Pipeline for Processing
def submit_pipeline (wdl_file, json_file):
	command_line = ('curl -s -X POST "%s"'
		' -H "accept: application/json" -H "Content-Type: multipart/form-data"'
		' -F "workflowSource=@%s"'
		' -F "workflowInputs=@%s;type=application/json"'
		% (cromwell_base, wdl_file, json_file))
	parsed_json = run_subprocess(command_line)
	cromwell_id = parsed_json['id']
	print ("Cromwell ID:  %s" % cromwell_id)
	return cromwell_id

# Check Cromwell Status of Specified ID
def check_status (cromwell_id):
	command_line = ('curl -s -X GET "%s/%s/status" -H "accept: application/json"'
		% (cromwell_base, cromwell_id))
	parsed_json = run_subprocess(command_line)
	status = parsed_json['status']
	return status

# Get Cromwell Pipeline VCF Output
def get_vcf_output(cromwell_id):
	command_line = ('curl -s -X GET "%s/%s/outputs"'
		' -H "accept: application/json"' % (cromwell_base, cromwell_id))
	parsed_json = run_subprocess(command_line)
	outputs = parsed_json['outputs']
	vcf_file = outputs["ecoliWorkflow.callVariants2.vcf_file"]
	return vcf_file

# Validate Results in the VCF File
def validate_vcf(vcf_file):
	fd = open (vcf_file)
	for line in fd:
		if not line.startswith("#"):
			variant_line = line
	fields = variant_line.split("\t")
	assert len(fields) > 2
	assert fields[0] == "gi|110640213|ref|NC_008253.1|"
	assert fields[1] == "736"
	assert fields[3] == "T"
	assert fields[4] == "G,C"

# Test the Pipeline (local, no docker)
def test_pipeline_local_no_docker():
	cromwell_id = submit_pipeline("ecoli_no_docker.wdl", "ecoli_local.json")
	status = "Running"
	while status == "Running" or status == "Submitted":
		time.sleep(5)
		status = check_status(cromwell_id)
		print (status)
	vcf_file = get_vcf_output(cromwell_id)[0]
	validate_vcf(vcf_file)