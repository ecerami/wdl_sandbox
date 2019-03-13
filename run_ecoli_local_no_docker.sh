python strip_docker.py ecoli.wdl > ecoli_no_docker.wdl
python create_inputs.py ecoli_base.json local > ecoli_local.json
java -jar womtool-38.jar validate ecoli_no_docker.wdl
java -jar cromwell-38.jar run ecoli_no_docker.wdl -i ecoli_local.json
