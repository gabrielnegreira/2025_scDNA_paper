#!/usr/bin/env bash
set -euo pipefail

#script to install ScisTree2 according to the steps described in https://github.com/yufengwudcs/ScisTree2/

#set a directory to install it
mkdir -p $HOME/opt
cd $HOME/opt

#clone the repo
if [ ! -d ScisTree2 ]; then
  git clone https://github.com/yufengwudcs/ScisTree2.git
fi
cd ScisTree2

#create a virtual environment for it
#get the base for conda
conda_base=$(conda info --base)
conda_sh="$conda_base/etc/profile.d/conda.sh"

#Verify it exists, then source it
if [[ -f "$conda_sh" ]]; then
  source "$conda_sh"
else
  echo "ERROR: Cannot find conda.sh at $conda_sh" >&2
  exit 1
fi

conda create --yes --name scistree2 python=3.9
conda activate scistree2

#install it
pip install --upgrade pip
pip install .
