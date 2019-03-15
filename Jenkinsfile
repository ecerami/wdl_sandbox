pipeline {
    agent { docker { image 'python:2.7.10' } }
    stages {
        stage('test') {
            steps {
                sh 'test.sh'
            }
        }
    }
}
