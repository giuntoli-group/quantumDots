# quantumDots
This repository contains the general LAMMPS codes used in our paper to simulate the assembly of monolayer quantum PbSe dot arrays. 


quantumDots-assembly.in is the general script used for all systems. Pressure and run time must be changed accordingly depending on the simulation.
Hexagonal.py is the python script used to generate the initial structures of 2D arrays with varying spacing orientation of the dots.
Single_QD contains the atomic coordinates for one individual quantum dot, used in the Hexagonal.py script to generate the entire structure.
