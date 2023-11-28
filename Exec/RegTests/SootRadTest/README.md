## SootRadTest
Testing of the coupling between PeleLMeX and PeleMP soot and radiation modules.

For now, PeleRad must be clone separately for this test case. In a convenient
location, run:

git clone https://github.com/AMReX-Combustion/PeleRad.git
export PELERAD_HOME=$(pwd)/PeleRad

Please specify the radiation database path in the input file before running.

echo "pelerad.kppath = "$PELERAD_HOME/data/kpDB/"" >> first-input
