import context

from py_cromwell.server import CromwellServer

# Test the E-Coli Pipeline
def test_ecoli_local_pipepline():
	server = CromwellServer()
	cromwell_id = server.submit_pipeline("ecoli_no_docker.wdl", "ecoli_local.json")
	vcf_file = server.get_output(cromwell_id, "ecoliWorkflow.callVariants2.vcf_file")
	validate_vcf(vcf_file)

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