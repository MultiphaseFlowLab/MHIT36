module load gcc
rm -r *.mod
make clean
rm -r output
mkdir output
make
./read_paraview
