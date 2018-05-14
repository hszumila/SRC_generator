# SRC_generator

This generator was made to look at kinematic distributions, especially for SRC in CLAS12. This code outputs a root file for checking the generated information as well as a LUND file (input.dat) which can be read into GEMC. The LUND file contains the incident and scattered electron, the proton and the recoil particle (default is assumed to be proton). 
```
 .x generatorSRC.C (beamP, zVtx, nEvents)
```
The beamP variable is the input beam momentum. The zVtx is the location of the generated event in the z direction where +z is downstream. The nEvents is the number of events you wish to generate. The generator assumes scattering on a carbon nucleus at rest. 