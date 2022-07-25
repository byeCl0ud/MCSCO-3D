# MCSCO-3D

## Overview
Small CLI program written in FORTRAN which works as a Monte Carlo spin crossover simulation 
in Ising model aproximation in three dimensions. 
It generates output file with 3 columns - temperature, energy and magnetization. 
Programs like gnuplot can be used to make graphs. 

![](https://github.com/byeCl0ud/MCSCO-3D/blob/development/MCSCO_animation3D.gif)

## Building
Use build script compile.sh to compile with Intel Fortran Compiler, 
or gfortran_compile.sh to compile with GNU Fortran: 

```sh
chmod +x compile.sh && ./compile.sh
chmod +x gfortran_compile.sh && ./gfortran_compile.sh
```

## Running
First set your input according to the input structure below in file in.2d. 
Then simply run compiled file: 
```sh
./a.out
```
## License
Published under the GNU GPL License. See `LICENSE.txt` for more information.
