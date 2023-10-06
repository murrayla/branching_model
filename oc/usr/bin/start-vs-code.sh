#!/bin/bash

# To access older or newer releases, go to https://code.visualstudio.com/updates/v1_51 and click the version of interest from the UPDATES list on the left-hand-side of the page and copy and paste the linux tarball link found at the top of the resulting page into the VSCODE_SOURCE variable below.
export VSCODE_VERSION=1.51.1
export VSCODE_SOURCE=https://update.code.visualstudio.com/${VSCODE_VERSION}/linux-x64/stable

# Specify the name of the folder that is extracted from the tarball.
export VSCODE_DIR=~/work/VSCode-linux-x64
# The portable version of VSCODE stores configuration files etc in a folder called data within the folder where the VSCODE tarball is extracted (see https://code.visualstudio.com/docs/editor/portable for more information).
export VSCODE_DATA_DIR=${VSCODE_DIR}/data

# Check if the Visual Studio executable is available.
FILE=${VSCODE_DIR}/code
if [ -f "$FILE" ]; then
    echo "Visual Studio Code is already installed."
else 
    echo "Installing Visual Studio Code."
    mkdir -p ${VSCODE_DIR}
    wget --output-document=${VSCODE_DIR}/installer.tgz ${VSCODE_SOURCE} \
    && cd ${VSCODE_DIR} || exit \
    && tar --strip-components=1 -xzf installer.tgz \
    && rm installer.tgz
    # Create the Visual Studio Code data directory
    mkdir ${VSCODE_DATA_DIR}
fi

# Execute Visual Studio Code (use --no-sandbox argument to workaround the following bug https://github.com/microsoft/vscode/issues/81358)
${VSCODE_DIR}/code --no-sandbox --no-xshm
