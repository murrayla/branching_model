#!/bin/bash

# Start OpenCMISS-iron virtual environment
. /home/jovyan/work/bashrc

export OpenCMISSLibs_DIR=/home/jovyan/opencmiss-build/opencmiss/install

export TESTING_DIR=~/work/functional-tests
mkdir -p ${TESTING_DIR}
cd ${TESTING_DIR} 
git clone https://github.com/OpenCMISS/functional_test_framework.git
git clone https://github.com/OpenCMISS/functional_test_database.git
mkdir functional_test_framework-build
cd functional_test_framework-build
cmake -DOpenCMISSLibs_DIR=${OpenCMISSLibs_DIR} -DTEST_DB=../functional_test_database/tests ../functional_test_framework
make

