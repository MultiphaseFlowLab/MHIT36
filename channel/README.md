# Channel configuratuion

Code for turbulent channel flow simulations (closed channel).
Solution of Navier-Stokes equations + phase-field method (ACDI) + energy equation (passive scalar).
For time integration, a fully explicit scheme is used: Euler explicit for ACDI + RK3 for Navier-Stokes and temperature

## Check list of features implemented
- Boundary condition for no-slip at the two walls ✅
- Laminar solution (uniform and non-uniform grids) ✅
- TDMA validation ✅ 
- Turbulent channel flow ✅ 
- Stretched grids ✅

**Autotuning of decomposition**  
- Default: `pr=0` and `pc=0` → enables **automatic process decomposition**.  
- Only input needed: **total number of MPI tasks**.  
- ✅ No recompilation required when changing MPI processes.  

**Optional features (conditional compilation)**  
- Phase-field module: Can be enabled or disabled. Default is single-phase.  
- Passive scaler: Can be enabled or disabled. Default is single-phase.  

## Turbulent channel flow 
- TCF at Rey_tau=590, Grid: 1536 x 768 x 576

![Test](../val/tcf2.png)


## Nodes numbering and staggered grid

X-pencil configuration:

![Test](../val/grid5.png)

## Contributing

We welcome all contributions that can enhance TCF36, including bug fixes, performance improvements, and new features. 
If you would like to contribute, please contact me or open an Issue in the repository.