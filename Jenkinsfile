node("conan") {
    stage("Checking out sources") {
        git branch: "jenkins-test", url: "https://github.com/dariusarnold/doublebeam/"
    }
    stage("Install dependencies with Conan") {
        sh "mkdir -p build && cd build && conan install .."
    }
    stage("Build") {
        cmakeBuild buildDir: 'build', buildType: 'Debug', generator: 'Ninja', installation: 'InSearchPath', steps: [[args: '--target doublebeam', withCmake: true]]
    }
}