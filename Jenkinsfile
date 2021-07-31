node("conan") {
    stage("Install dependencies with Conan") {
        echo "Installing dependencies..."
        sh "mkdir -p build && cd build && conan install .."
    }
    stage("Build") {
        echo "Building..."
    }
}