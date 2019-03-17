echo Creating Input JSON File
python bin/create_inputs.py wdl/ecoli_base.json google > wdl/ecoli_google.json
echo Validating WDL
java -jar womtool-38.jar validate wdl/ecoli.wdl
echo Running Pipeline via Cromwell + Google Genomics API
java -Dconfig.file=conf/google.conf -jar cromwell-38.jar run wdl/ecoli.wdl -i wdl/ecoli_google.json
