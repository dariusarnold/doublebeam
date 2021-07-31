pipeline {
    agent {label "conan"}
    stages {
        stage("Checking out sources") {
            steps {
                git branch: "jenkins-test", url: "https://github.com/dariusarnold/doublebeam/"
            }
        }
        stage("Install dependencies with Conan") {
            steps {
                echo "Installing dependencies..."
                sh "mkdir -p build && cd build && conan install .."
            }
        }
        stage("Build") {
            steps {
                echo "Building..."
            }
        }
    }
}