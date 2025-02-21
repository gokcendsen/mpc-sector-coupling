
# mpc-sector-coupling
Model Predictive Control (MPC) for cost-optimal predictive operation of sector-coupled energy plants. Developed for our SENSYS 2025 paper.

# **SENSYS Energy Management Code**  

This repository contains the Julia implementation of a **Model Predictive Control (MPC)-based energy management system** for sector-coupled energy plants.  

## 📂 **Files**  
- **`SENSYS_MPC_optimization.jl`** – Julia script implementing the MPC optimization.  
- **`SENSYS_energy_data.xlsx`** – Excel file containing energy demand data used in simulations.  

## 🚀 **Usage**  
Run the Julia script to execute the optimization:  
```julia
julia SENSYS_MPC_optimization.jl
```

📌 Important  
Ensure that `SENSYS_energy_data.xlsx` is in the same directory as `SENSYS_MPC_optimization.jl` before running the script.

📌 Reference
This code is associated with the paper:
"A cost-optimal predictive operation scheme for sector-coupled energy plants with start-up delays and start-up costs" (submitted to SENSYS 2025).

📝 License
This project is licensed under the MIT License – see the LICENSE file for details.

📖 Cite this work
[![DOI](https://zenodo.org/badge/936695872.svg)](https://doi.org/10.5281/zenodo.14906886)









