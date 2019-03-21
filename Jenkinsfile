pipeline {
  agent {
    docker {
      image 'wolfram-docker-11.3.0'
    }

  }
  stages {
    stage('Run tests') {
      steps {
        dir(path: 'ReggeWheeler') {
          checkout scm
          sh 'Tests/AllTests.wls'
          junit 'TestReport.xml'
        }

      }
    }
  }
  options {
    skipDefaultCheckout(true)
  }
}