node("conan") {
    stage("Checking out sources") {
        checkout poll: false, scm: [$class: 'GitSCM', branches: [[name: '*/jenkins-test']],
        browser: [$class: 'GithubWeb', repoUrl: 'https://github.com/dariusarnold/doublebeam/'],
        extensions: [[$class: 'CleanCheckout']], userRemoteConfigs: [[url: 'https://github.com/dariusarnold/doublebeam.git']]]
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