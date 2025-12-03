# Channel configuratuion

Code for turbulent channel flow simulations (closed channel).
Solution of Navier-Stokes equations + phase-field method (ACDI) + energy equation (passive scalar).
For time integration, a fully explicit scheme is used: Euler explicit for ACDI + RK3 for Navier-Stokes and temperature

## Main features at a glance:
- No-slip boundary conditions at the two walls ($z=\pm h$)
- Poisson solver based on FFT2 along $x$ and $y$ + TDMA along $z$.
- Non-uniform grid along $z$, any stretching function can be used (defined in read_input.f90)
- Domain of arbitary dimensions along all three directions

**Autotuning of decomposition**  
- Default: `pr=0` and `pc=0` → enables **automatic process decomposition**.  
- Only input needed: **total number of MPI tasks**.  
- No recompilation required when changing MPI processes.  

**Optional features (conditional compilation)**  
- Phase-field module: Can be enabled or disabled. Default is single-phase.  
- Passive scaler: Can be enabled or disabled. Default is single-phase.  

## Turbulent channel flow 
- Turbulent channel flow at $Re_\tau=590$, Grid: $1536 \times 768 \times 576$

![Test](../val/tcf2.png)

## Nodes numbering and staggered grid

X-pencil configuration:

![Test](../val/grid_channel.png)
