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
    stage("Build tests") {
        cmakeBuild buildDir: 'build', buildType: 'Debug', generator: 'Ninja', installation: 'InSearchPath', steps: [[args: '--target Unit_Tests_run', withCmake: true]]
    }
    stage("Run tests") {
        dir("build") {
            sh "pwd && ./tests/unit_tests/Unit_Tests_run"
        }
    }
}

node("conan") {
    stage("Clean workspace") {
        echo "Deleting build artefacts"
        //cleanWs deleteDirs: true
    }
}