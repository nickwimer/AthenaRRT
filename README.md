# AthenaRRT
Athena 4.2 with support for reacting Rayleigh Taylor instability.

Boundary Conditions:  Currently the sponge boundary conditions are hard coded in integrate_2d.c so they do not need to be activated in the athinput file.  However there are several variables that need to be in athinput for the boundary conditions to function.

ASZx  --> the size of the sponge domain in the x direction (in physical dimensions)
ASZy  --> the size of the sponge domain in the y direction
sig0x --> damping coefficient of the sponge layer in the x direction (should be 0 if using periodic boundaries in the x direction)
sig0y --> damping coefficient of the sponge layer in the y direction (try 50 to start)
