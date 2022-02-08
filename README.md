# TeraHertz radiation from laser pumped heterostructures

This repository's purpose is to model and simulate laser pumped heterostructures. For detailed description of the modelling, 
please check the EFiieldSimulation.Simulation module, or the associated report.

### Packages

##### Fortran Simulation

The 

##### EFieldSimulation
Represent the simulated system. It takes the output from the Fortran code as its input, and outputs E-fields.

This module's purpose is to simplify the execution of Pavel Balaz's Superdiffusive Spin Transport script in Fortran. 

EfieldSimulation:
    

Frequency:
    This module contains script that analyze the frequenzy properties of the simulated E-field.

PipelineTools:
    This module contains pipeline functions that helps perform the entire Python simulation, from spin current to final results
     