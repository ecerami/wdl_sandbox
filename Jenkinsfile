pipeline {
    agent { docker { 
            image 'ecerami/wdl_test:latest'
            args '-v /Users/ecerami/dev/wdl_sandbox/cromwell-executions/ecoliWorkflow/:/Users/ecerami/dev/wdl_sandbox/cromwell-executions/ecoliWorkflow/'
    } }
    stages {
        stage('test') {
            steps {
                sh 'bin/test.sh'
            }
        }
    }
    post {
        always {
            junit 'junit/*.xml'
        }
    }
}
