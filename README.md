---
title: CFS-PML2nd Manual
date: 2020-10
author: Yufeng Wang
mathjax: true
---


## Overview of `CFS-PML2nd` package

`CFS-PML2nd` is a C-based code package that implements unsplit CFS-PML for fractional viscoacoustic simulation. This package is provided for verifying the feasibility and stability of the proposed CFS-PML when used to the second-order rather than first-order wave equation. We provide two package versions: `cfspml2nd` for viscoacoustic modeling with CFS-PML; `pml2nd` for viscoacoustic modeling with standard PML. The package has been well test on Linux OS (Ubuntu 18.04 TLS) with `g++` and `SeismicUnix` (<https://github.com/JohnWStockwellJr/SeisUnix>) available.

## The architecture of `CFS-PML2nd` package 

-   `input`: accurate velocity and $Q$ model for $Q$-RTM:
    - `acc_vp.dat`: quasi-Marmousi velocity model;
    - `acc_Qp.dat`: quasi-Marmousi $Q$ model;
-   `output`: generated results include:(seismograms, snapshots, and energydecay). 
    - `visco_1_1_profile.dat`: seismograms recorded at the surface;
    - `energydecay.dat`: energy decay of the computational domain;
    - `visco_1_1_receive(1-4).dat`: seismic records at four receivers;
    - `visco_1_1_snapshots1_(0-8000).dat`: snapshots recorded at every 100 time steps;
-   `plot`: scripts for plotting figures, which includes:
    - `energydeacy.m`: plot the energy decay of the CFS-PML and standard PML schemes;
    - `snapscpml.m`: plot the snapshots of the CFS-PML scheme at some time slices;
    - `snapspml.m`: plot the snapshots of the standard PML scheme at some time slices.
-   `Myfunctions.h`: header file;
-   `Myfunctions.c`: subfuncations file;
-   `viscocpml.c`: C code file, there are serveal important parameters to control performance of `CFS-PML2nd`, which includes:
    - `nt`: the time steps for simulation, if you want to have a quick test, set `nt=2000` or less, 
``` c
	int nt=8000;
```
    - `Memory_len` you can set `Memory_len=40` for quick test, and `Memory_len=300` or more for removing the side effects of this approximation.
``` c
int Memory_len=300;	//>=4
```
    - `Gmax` and `Amax`: these parameters are defined for CFS-PML, one can change them for test. If `Gmax = 1.0` and `Amax = 1.0*pi*f0`, the SFS-PML reduces to standard PML scheme.
``` c
	float R=1E-6;
	float Vmax=3000.0;
	float Gmax=1.0;
	float Amax=1.0*pi*f0;
```
   -`Makefile`: excution script.

## Prerequisites

`CFS-PML2nd` package is developed under `Linux` system, which should be equipped with the following environments:

- `C` environment (`gcc` or `g++`);
- `SeismicUnix`;
- `Matlab`;


## How to run this package

- Step 1: Confirm the environment in `Makefile`, and replace the folder path with your own enviroment path; 

``` bash
INC1=/home/wyf/SeismicUnix/include
INC2=/home/wyf/SeismicUnix/src/Complex/include
LIK1=/home/wyf/SeismicUnix/lib
LIK2=/home/wyf/SeismicUnix/src/Complex/lib
LIB=-lsu -lpar -lcwp -lm

all: viscocpml.c
	g++ -o a.out viscocpml.c -fopenmp -I$(INC1) -I$(INC2) -L$(LIK1) -L$(LIK2)  $(LIB)
	clear        
	nohup ./a.out&
```

- Step 2: Run the `Makefile` by the command line: `make`;
- Step 3: View generated files in the folder `./ouput`;
- Step 4: Plot figures by run `/plot/energydeacy.m`, `/plot/snapscpml.m`, `/plot/snapspml.m`.

## Contact me

I am Yufeng Wang, an associate researcher from China University of Geosciences, Wuhan. If you have any question about this coda package, please feel free to contact me by [Email:wangyufeng@cug.edu.cn](wangyufeng@cug.edu.cn).

## Copyright

`CFS-PML2nd` is a C-based code package that implements unsplit CFS-PML scheme for fractional viscoacoustic simulation.

The Marmousi model used in this package is available for download from `Madagascar` MainPage <http://www.ahay.org/data/marm2>

Copyright (C) 2020  China University of Geosciences, Wuhan (Yufeng Wang)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
