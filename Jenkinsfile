node("conan") {
    stage("Checking out sources") {
        checkout git branch: "jenkins-test"
    }
    stage("Install dependencies with Conan") {
        echo "Installing dependencies..."
        sh "mkdir -p build && cd build && conan install .."
    }
    stage("Build") {
        echo "Building..."
    }
}