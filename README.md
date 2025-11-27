# MHIT36

Multi-GPU code for intreface-resolved simulations of multiphase turbulence in homogeneous isotropic turbulence and channel flow configurations.
The code relies on direct numerical simulation of Navier-Stokes equations coupled with a phase-field method (ACDI) for interface description.
Tracking of Lagrangian particles (tracers) is also supported.
The code parallerelization relies on the cuDecomp library.


~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████       
████  ████ ██   ██ ██    ██         ██ ██              
██ ████ ██ ███████ ██    ██     █████  ███████   
██  ██  ██ ██   ██ ██    ██         ██ ██    ██     
██      ██ ██   ██ ██    ██    ██████   ██████        
~~~


If you use this code, please cite the following works: 
```bibtex
  @article{roccon2025,
  title   = {MHIT36: A Phase-Field Code for Gpu Simulations of Multiphase Homogeneous Isotropic Turbulence},
  author  = {Roccon, A. and Enzenberger, L. and Zaza, D. and Soldati, A.},
  journal = {Computer Physics Communications},
  year    = {2025},
  volume  = {314},
  issue   = {109804},
  doi     = {https://doi.org/10.1016/j.cpc.2025.109804}
}
```

```bibtex
  @article{roccon2026,
  title   = {MHIT36: Extension to wall-bounded turbulence and scalar transport equation},
  author  = {Roccon, A.},
  journal = {Computer Physics Communications},
  year    = {2026},
  volume  = {---},
  issue   = {109956},
  doi     = {https://doi.org/10.1016/j.cpc.2025.109956}
}
```

![Test](val/render2.jpg)


## How to run the code
- **Compile the cuDecomp library** using `*_lib.sh`.  
  - The resulting modules and libraries will be located in:  
    - `cuDecomp/build/lib`  
    - `cuDecomp/build/include`  
  - Alternatively, you can compile manually using the **CMake** setup provided in the cuDecomp repository.  
- **Verify the cuDecomp build** to ensure it compiled correctly (must be built using **HPC-SDK**).  
- **Folder `hit/` or `channel/`**: contains the source code of the multi-GPU version of the code for homogeneous isotropic and channel flow configurations, respectively. 
  - Use one of the scripts (`local.sh`, `leo.sh`, or `mn5.sh`) to compile and run.  
- **Autotuning (default mode):**  
  - By default, the code is set to full autotuning (`pr=0`, `pc=0`).  
  - cuDecomp will automatically determine the best decomposition at runtime, based only on the total number of tasks.  
  - This means recompilation is **not required** when changing the number of MPI processes.  
  - The backend is also selected automatically.  
- **Phase-field module (optional):**  
  - Controlled via a conditional compilation flag.  
  - By default, only the single-phase version is enabled.  
- **Scalar equation module (optional):**  
  - Controlled via a conditional compilation flag.  
  - By default, only the single-phase version is enabled. 
- **Particles module (optional):**  
  - Controlled via a conditional compilation flag.  
  - This part is supported only for the hit configurations ATM.
- **Additional information**
  - Additional information on the numerical implementatioon, nodes collocation and validation are available in the respective README.md file present in [/hit](hit) and [/channel](channel) folders. 


## Reference performance and scaling (hit configuration)
Performance (NS only)
* 128 x 128 x 128    |   1 x A100@Leonardo  |   1 ms/timestep 
* 256 x 256 x 256    |   1 x A100@Leonardo  |   8 ms/timestep
* 512 x 512 x 512    |   1 x A100@Leonardo  |  65 ms/timestep 
* 128 x 128 x 128    |   4 x A100@Leonardo  |   1 ms/timestep
* 256 x 256 x 256    |   4 x A100@Leonardo  |   3 ms/timestep
* 512 x 512 x 512    |   4 x A100@Leonardo  |  18 ms/timestep 
* 512 x 512 x 512    |   4 x H100@MN5-ACC   |  14 ms/timestep 
* 1024 x 1024 x 1024 |   4 x A100@Leonardo  | 150 ms/timestep 
* 2048 x 2048 x 2048 |  64 x A100@Leonardo  | 330 ms/timestep
* 4096 x 4096 x 4096 | 256 x A100@Leonardo  | 780 ms/timestep

## Reference performance and scaling (channel configuration)
Performance (NS only)
* 256 x 128 x 200    |   2 x RTX5000 16GB   |  31 ms/timestep 
* 512 x 256 x 384    |   4 x A100 64 GB     |  15 ms/timestep
* 1536 x 768 x 576   |   4 x A100 64 GB.    | 220 ms/timestep
* 2048 x 768 x 576   |   4 x A100 64 GB     | 323 ms/timestep
* 2048 x 768 x 576   |  16 x A100 64 GB     | 127 ms/timestep
* 3456 x 1296 x 960  |  16 x A100 64 GB     | 720 ms/timestep

- Phase-field introduces about 20% of overhead compared to NS only
- Max resolution tested: 4096^3
- Max number of GPUs used: 512 (Leonardo) and 1024 (MN5)


## Validation

Benchamrk present in [VanRees2021](https://www.sciencedirect.com/science/article/abs/pii/S0021999110006467) and [Costa2018](https://www.sciencedirect.com/science/article/pii/S089812211830405X).

Time evolution of the viscous dissipation:

![Test](val/val.png)

Single refers to the single-GPU implmenetation of the code, which is available [here](https://github.com/aroccon/MHIT36_single) and can be used as a sandbox for testing being very easy to modify. Other validations tests are present in the two CPC papers for both hit and channel configurations. See also the supplementary material.

## Contributing

We welcome all contributions that can enhance MHIT36, including bug fixes, performance improvements, and new features. 
If you would like to contribute, please contact aroccon or open an Issue in the repository.

## Acknowledgements

We would like to thank the follwing people for the support received during the code development and optimization as well as CINECA and the Open Hackathons program:
- Matt Bettencourt
- Josh Romero
- Laura Bellentani
- Alessio Piccolo
- Pedro Costa 
- Simone Di Giorgio