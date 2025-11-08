make clean
rm -r *.mod
make &> /dev/null
make
./stats
