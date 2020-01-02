#!/bin/bash
systemName=$(uname -a)

echo =======================================================================================
echo Install Solver for Climate Project
echo =======================================================================================

echo Starting installation process...
if [[ $systemName == *"Darwin"* ]]; then
  export MACOSX_DEPLOYMENT_TARGET=10.9
fi
echo ===============================================================================
echo Step 1: Install numba and pybind11
pip install pybind11
pip install plotly

echo ===============================================================================
echo Step 2: Install model solution core
pip install ./src/cppcore
