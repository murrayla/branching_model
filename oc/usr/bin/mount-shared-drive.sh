#!/bin/bash

sudo mkdir -p /shared_drive
sudo mount -v -t cifs //IP_ADDRESS/DRIVE_PATH/ /shared_drive -o username=USER,domain=uoa,vers=3.0
