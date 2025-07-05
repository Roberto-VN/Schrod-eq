# Schrödinger equation

## Introduction
In this study we have made a numerical approach to the two dimensional time-dependent Schrödinger equation in order to study the wavy behaviour of a particle governed by its wave function. Specifically, in a double-slit setup and some variations.

## Functions
The files `Functions.hpp` & `Functions.cpp` define the functions that will be requiered for the main code. These functions will create the initial state of the wave, and will be able to advance the wave along the time for different settings.  (For more specific infomation about the functions, see the code).

### Note
The parameters of the simulation are defined in the terminal. These are the inputs: h (distance variation), Δt (time variation), T (duration of the simulation), x<sub>c</sub> (x-position of the coordinates of the centre of the initial wave packet), σ<sub>x</sub> (x initial width of the wave packet), p<sub>x</sub> (x wave packet momenta), y<sub>c</sub> (y-position of the coordinates of the centre of the initial wave packet), σ<sub>y</sub> (y initial width of the wave packet), p<sub>y</sub> (x wave packet momenta), v<sub>0</sub> (potential of the wall), slits (number of slits), mode. The last parameter defines the mode which can takes values from 1 to 4. Depending on this value, the code will return different information of our system. 

## Probability deviation
The particle is supose to have a probability of 1 of being found inside the box. However, because of the numerical aproximations there will be a deviation of this value. The following commands with mode = 1 return this deviation along the time in a file (Prov_dev.txt). 

The following build command just need to be use one time and will work for every run command on this README.

### Build:
```bash
g++ -o Schrodinger Schrodinger.cpp Functions.cpp -larmadillo -llapack -lblas
```

### No barrier case

```bash
./Schrodinger 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.1 0 0 2 1
```

### Double slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.1 0 10e10 2 1
```

## Specific wave values 
The mode = 2 returns: The probabilities of finding the particle in each point at t = {0, T/2, T} (File: Prob_map.txt); values of the wave in each point at t = {0, T/2, T} (File: Wave_val.txt); probabilities of findin the particle along the y axis when the particle has been meassured in a certain x  at t = T (File: Prob_y.txt). It is important to know that the x position where the particle is meassure must be changed in the code, not in the terminal.

### Single slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.2 0 10e10 1 2
```

### Double slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.2 0 10e10 2 2
```

### Triple slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.2 0 10e10 3 2
```

## Simulation 
The mode = 3 returns a file (Wave_sim.txt) with the probabilities of finding the particle in each point of the box at each time step.

### Single slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.2 0 10e10 1 3
```

### Double slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.2 0 10e10 2 3
```

### Triple slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.2 0 10e10 3 3
```

## Wave values
With mode = 4, the output will be a file (Schr_eq.txt) with the values of the wave in each point of the box at each time step.

### Single slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.2 0 10e10 1 4
```

### Double slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.2 0 10e10 2 4
```

### Triple slit case

### Run:
```bash
./Schrodinger 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.2 0 10e10 3 4
```
