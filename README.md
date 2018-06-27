#About ItrAs

This is a software tool relying on GATB-CORE library. It essentially copies the theory behind IDBA-UD

It uses the backbone of Minia though to do iterative kmer assemly. You specify and minimum k, a maximum k, and a step size. The assembler will then build a graph for the mink size and produce contigs. The contigs and the read file are then used to build the next graph where new contigs are produced. These contigs along with reads are used to build the next graph... and so on. 

It has Karect Read Error correction integrated (enabled with -pre_correct 1). Karect is limited to the free amount of RAM on the machine you run on.

To build, just type 'make'

To Test, type 'make test' - This will run the program on the GAGE Staphylococcus aureus AllPaths corrected dataset. (located in tests folder)

To run, specify the minimum k, maximum k, step size, input file and output directory. 
./build/tools/ItrAs --help (to show all options)

Example:
./build/tools/ItrAs -in inputfile.fa -mink 35 -step 10 -maxk 95 -out OutputDir
