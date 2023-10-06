#!/bin/bash
export PYCHARM_BUILD=2020.2.3
export PYCHARM_SOURCE=https://download.jetbrains.com/python/pycharm-community-${PYCHARM_BUILD}.tar.gz

# Specify the name of the folder that is extracted from the tarball.
export PYCHARM_DIR=~/work/PyCharm

FILE=${PYCHARM_DIR}/bin/pycharm.sh
if [ -f "$FILE" ]; then
    echo "Pycharm already installed."
else 
    echo "Installing pycharm."
    mkdir -p ${PYCHARM_DIR}
    wget --output-document=${PYCHARM_DIR}/installer.tgz ${PYCHARM_SOURCE} \
    && cd ${PYCHARM_DIR} || exit \
    && tar --strip-components=1 -xzf installer.tgz \
    && rm installer.tgz
fi

~/work/PyCharm/bin/pycharm.sh


