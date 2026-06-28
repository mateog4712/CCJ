# CCJ

### Description
Contains the software implementation for CCJ-PF.
CCJ-PF is a partition function based method for predicting the psuedoknotted secondary structures of RNA sequences.      

CCJ-PF should work on most Linux or Mac machines.       

### Installation:  
Requirements: A compiler that supports C++17 standard (tested with g++ version 4.7.2 or higher)  and CMake version 3.1 or greater.    

[CMake](https://cmake.org/install/) version 3.1 or greater must be installed in a way that CCJ can find it.    
To test if your Mac or Linux system already has CMake, you can type into a terminal:      
```
cmake --version
```
If it does not print a cmake version greater than or equal to 3.1, you will have to install CMake depending on your operating system.

#### Mac:    
Easiest way is to install homebrew and use that to install CMake.    
To do so, run the following from a terminal to install homebrew:      
```  
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"   
```    
When that finishes, run the following from a terminal to install CMake.     
```   
brew install cmake   
``` 
#### Linux:    
Run from a terminal     
```
wget http://www.cmake.org/files/v3.8/cmake-3.8.2.tar.gz
tar xzf cmake-3.8.2.tar.gz
cd cmake-3.8.2
./configure
make
make install
```
[Linux instructions source](https://geeksww.com/tutorials/operating_systems/linux/installation/downloading_compiling_and_installing_cmake_on_linux.php)

#### Steps for installation   
1. [Download the repository](https://github.com/mateog4712/CCJ/archive/master.zip) and extract the files onto your system.
2. From a command line in the root directory (where this README.md is) run
```
cmake -H. -Bbuild
cmake --build build
```   
If you are getting errors about your comiler not having C++11 features, you may need to specify a specific compiler, such as g++.
If you want to do so you can instead run something like   
```
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
```   

Help
========================================

```
Usage: CCJ [options] [input sequence]
```

Read input file from cmdline; predict minimum free energy and optimum structure using the time- and space-efficient MFE RNA folding algorithm.

```
-h, --help                Print help and exit
-V, --version             Print version and exit
-i, --input-file=STRING   Give a path to an input file containing the sequence (and input structure if known)
-o, --output-file=STRING  Give a path to an output file which will the sequence, and its structure and energy
-d, --dangles=INT         Specify the dangle model to be used (default=`2')
-P, --paramFile=STRING    Read energy parameters from paramfile, instead of using the default parameter set.
-s, --samples=INT         Give the number of samples for the stochastic backtracking (default=`1000')
-f, --fatgraph=INT        Give the number of fatgraphs outputted, along with their frequencies (default=`1')
-O, --print-samples       Print the samples with their multiplicities (default=off)",
    --noConv              Do not convert DNA into RNA. This will use the Matthews 2004 parameters for DNA (default=off)
    --noGU                Turn off G-U and U-G (and G-T and T-G) base pairing (default=off)
    --noPS                Don't create a Postscript drawing of the base pair probabilities (default=off)
```

Remarks:
    The default parameter file is DP09. This can be changed via -P and specifying the parameter file you would like

#### Example:
    assume you are in the CCJ directory
    ./build/src/CCJ CGCGAGGGGCGCGAGGGGCCCCCCCC
    ./build/src/CCJ -P "params/rna_Turner04.par" CGCGAGGGGCGCGAGGGGCCCCCCCC
    ./build/src/CCJ -s 100000 -P "params/rna_Turner04.par"
    ./build/src/CCJ -i Examples/Cases.fa


### Licence    
CCJ-PF is copyrighted under GNU General Public Licence.

### Disclaimer
Although the authors have made every effort to ensure that CCJ-PF correctly implements the underlying models and fullfills the functions described in the documentation, neither the authors nor the University of Alberta guarantee its correctness, fitness for a particular purpose, or future availability.

### Contact  
If you have any issues or feature requests, please contact Mateo Gray: mateo2@ualberta.ca
