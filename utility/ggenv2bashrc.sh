#!/bin/bash

echo "# Use Graziano's environment" >> ~/.bashrc
echo "# Use Graziano's environment" >> ~/.tcshrc
echo "source /home/netapp-clima/users/ggiulian/minter-19.sh"   >> ~/.bashrc
echo "source /home/netapp-clima/users/ggiulian/minter-19.csh"  >> ~/.tcshrc

source ~/.bashrc
