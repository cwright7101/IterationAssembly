#About ItrAs

This is a software tool relying on GATB-CORE library. It essentially copies the theory behind IDBA-UD

It uses the backbone of Minia though to do iterative kmer assemly. 

To build, just type 'make'

To run, specify the minimum k, maximum k, step size, input file and output directory. 
./build/tools/ItrAs --help (to show all options)

Example:
./build/tools/ItrAs -in inputfile.fa -mink 35 -step 10 -maxk 95 -out OutputDir
