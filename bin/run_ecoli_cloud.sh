python create_inputs.py ecoli_base.json google > ecoli_google.json
java -jar womtool-38.jar validate ecoli.wdl
java -Dconfig.file=google.conf -jar cromwell-38.jar run ecoli.wdl -i ecoli_google.json
