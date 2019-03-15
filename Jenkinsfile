pipeline {
    agent { docker { image 'grihabor/pytest:python3.7-alpine' } }
    stages {
        stage('test') {
            steps {
                sh './test.sh'
            }
        }
    }
}
