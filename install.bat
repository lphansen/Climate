
echo off
echo =======================================================================================
echo Install Solver for Climate Project
echo Questions: please contact Macro Finance Research \(MFR\) Team at Becker Friedman Institute
echo URL: https://stackoverflow.com/questions/1098786/run-bash-script-from-windows-powershell
echo =======================================================================================

echo Starting installation process...

echo ===============================================================================
echo Step 1: Install plotly and pybind11
pip install pybind11
pip install plotly

echo ===============================================================================
echo Step 2: Install model solution core
pip install ./src/cppcore
