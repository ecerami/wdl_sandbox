echo Stripping out Docker from WDL
python bin/strip_docker.py wdl/ecoli.wdl > wdl/ecoli_no_docker.wdl
echo Creating Input JSON File
python bin/create_inputs.py wdl/ecoli_base.json local > wdl/ecoli_local.json
echo Validating WDL
java -jar womtool-38.jar validate wdl/ecoli_no_docker.wdl
echo Running Pipeline via Cromwell Local
java -jar cromwell-38.jar run wdl/ecoli_no_docker.wdl -i wdl/ecoli_local.json
