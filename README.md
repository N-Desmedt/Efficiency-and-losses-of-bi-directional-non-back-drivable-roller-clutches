# Efficiency and Losses of Bi-Directional, Non-Back-Drivable Roller Clutches

This repository contains the code and dataset for the paper:  
**“Efficiency and Losses of Bi-Directional, Non-Back-Drivable Roller Clutches.”**

---

## Repository Contents

### `Loss_computation.m`
A MATLAB script for computing the losses across all configurations of the clutch and the transitions between them.

### `Dataset.mat`
A MATLAB dataset containing the following variables:

---

### **Load Experiment Variables**

- **`Load_input_torque`**  
  Torque (Nm) measured on the input sensor during the load experiment, after filtering and only in steady state.  
  Rows correspond to each speed: *5 rpm, 10 rpm, 15 rpm, 20 rpm, 30 rpm, 40 rpm, 50 rpm, 60 rpm*.  
  Columns correspond to each data point collected.

- **`Load_output_torque`**  
  Torque (Nm) measured on the output sensor during the load experiment, after filtering and only in steady state.  
  Rows correspond to each speed: *5 rpm, 10 rpm, 15 rpm, 20 rpm, 30 rpm, 40 rpm, 50 rpm, 60 rpm*.  
  Columns correspond to each data point collected.

- **`Load_speed`**  
  Speed (rad/s) during the load experiment, after filtering and only in steady state.  
  Rows correspond to each experiment at the different speeds listed above.  
  Columns correspond to each data point collected.

- **`Load_mean_speed`**  
  Measured mean speed (rad/s) during each load experiment.

- **`Load_mean_loss`**  
  Mean loss (W) during the load experiment at each speed.

- **`Load_std_loss`**  
  Standard deviation of the loss during the load experiment at each speed.

---

### **No-Load Experiment Variables**

- **`No_load_input_torque`**  
  Torque (Nm) measured on the input sensor during the no-load experiment, after filtering and only in steady state.  
  Rows correspond to each speed: *0 rpm, 100 rpm, 200 rpm, 300 rpm, 400 rpm, 500 rpm, 750 rpm, 1000 rpm, 1250 rpm, 1500 rpm, 1750 rpm, 2000 rpm*.  
  Columns correspond to each data point collected.

- **`No_load_speed`**  
  Speed (rad/s) during the no-load experiment, after filtering and only in steady state.  
  Rows correspond to each experiment at the speeds listed above.  
  Columns correspond to each data point collected.

- **`No_load_mean_speed`**  
  Measured mean speed (rad/s) during each no-load experiment.

- **`No_load_mean_loss`**  
  Mean loss (W) during the no-load experiment at each speed.

- **`No_load_std_loss`**  
  Standard deviation of the loss during the no-load experiment at each speed.

---

## Citation
If you use this dataset or code, please cite the corresponding paper: “Efficiency and Losses of Bi-Directional, Non-Back-Drivable Roller Clutches.”
