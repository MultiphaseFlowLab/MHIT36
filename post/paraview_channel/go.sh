make clean
rm -r *.mod
make &> /dev/null
rm -r output
mkdir output
make
./read_paraview
