pipeline {
    agent { docker { 
            image 'ecerami/wdl_test:latest'
            args '-v /Users/ecerami/dev/wdl_sandbox/cromwell-executions/ecoliWorkflow/:/Users/ecerami/dev/wdl_sandbox/cromwell-executions/ecoliWorkflow/4'
    } }
    stages {
        stage('test') {
            steps {
                sh './test.sh'
            }
        }
    }
}
