import subprocess
import shlex
import json
import time
import os

# Cromwell Server -- Python Client
class CromwellServer:
	def __init__(self):
		if os.environ['JENKINS_WDL'] == "true":
			self.__cromwell_base = "http://docker.for.mac.localhost:8000/api/workflows/v1"
		else:
			self.__cromwell_base = "http://localhost:8000/api/workflows/v1"

	# Run a subprocess and capture the JSON Response
	def __run_subprocess(self, command_line):
		args = shlex.split(command_line)
		sub = subprocess.check_output(args)
		parsed_json = json.loads(sub)
		return (parsed_json)

	# Submit new WDL Pipeline for Processing
	def submit_pipeline (self, wdl_file, json_file):
		command_line = ('curl -s -X POST "%s"'
			' -H "accept: application/json" -H "Content-Type: multipart/form-data"'
			' -F "workflowSource=@%s"'
			' -F "workflowInputs=@%s;type=application/json"'
			% (self.__cromwell_base, wdl_file, json_file))
		parsed_json = self.__run_subprocess(command_line)
		cromwell_id = parsed_json['id']

		# Wait until completed
		status = "Running"
		while status == "Running" or status == "Submitted":
			time.sleep(5)
			status = self.check_status(cromwell_id)
		return cromwell_id

	# Check Cromwell Status of Specified Run ID
	def check_status (self, cromwell_id):
		command_line = ('curl -s -X GET "%s/%s/status" -H "accept: application/json"'
			% (self.__cromwell_base, cromwell_id))
		parsed_json = self.__run_subprocess(command_line)
		status = parsed_json['status']
		return status

	# Get Cromwell Output
	def get_output(self, cromwell_id, target_file):
		command_line = ('curl -s -X GET "%s/%s/outputs"'
			' -H "accept: application/json"' % (self.__cromwell_base, cromwell_id))
		parsed_json = self.__run_subprocess(command_line)
		outputs = parsed_json['outputs']
		vcf_file = outputs[target_file][0]
		return vcf_file