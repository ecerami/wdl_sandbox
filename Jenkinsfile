pipeline {
    agent { docker { image 'ecerami/wdl_test:latest' } }
    stages {
        stage('test') {
            steps {
                sh './test.sh'
            }
        }
    }
}
