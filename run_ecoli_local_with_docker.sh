python create_inputs.py ecoli_base.json local > ecoli_local.json
java -jar womtool-38.jar validate ecoli.wdl
java -jar cromwell-38.jar run ecoli.wdl -i ecoli_local.json
