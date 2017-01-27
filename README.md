# NetworkPercolation
Network Percolation Analysis

Hello!

This is a Fortran90 code to calculate cluster sizes using network percolation analysis in a molecular trajectory. This code uses molecular center of masses and periodic boundary conditions (input parameter) for the calculation and uses unformatted files (gromacs: xtc). The part that reads the xtc format was written by James W. Barnett jbarnet4@tulane.edu. The code has comments all over, but feel free to reach out to me at leejunta@grinnell.edu if you have any more questions. Special thanks to Roger Rousseau (PNNL) for assisting me in developing this code.

Here's how to run the code:

1) Install the xtc-library (http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library)

2) Copy the xtc-interface.f90 file into the appropriate folder

3) Edit the shell script to input the information you want inputted into the fortran code

4) gfortran -c xtc-interface.f90 -lxdrfile

5) gfortran xtc-interface.o network_percolation.f90 -lxdrfile -o network_percolation.x

6) chmod +x network_percolation.sh

7) ./network_percolation.sh

Source code available at https://github.com/leejunta/NetworkPercolation/
