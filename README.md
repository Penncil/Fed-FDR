# Documentation for "Federated Feature Selection with False Discovery Rate Control"

## Simulation Studies

###  Implementation of simulation under different settings:

#### Folder name: 
- simulation_result

#### R scripts to run:
 - simulation_n500p500.R
 - simulation_n500p1000.R
 - simulation_n1000p500.R
 - simulation_scalebility.R

#### Instruction for replication of simulation studies:
1. Download the R scripts from our repo
2. Use R or RStudio to run the R scripts
3. Results will be automatcially saved as .rds files
4. To generate the simulation results figures as in the manuscript, please download and run the following corresponding R scripts:
   
   4.1. Figure1.R
   
   4.2. Figure2.R
   
   4.3. Figure3.R


## Data Application Code

###  Running the code on real-world data

#### Folder name:
- use case

#### R script to run:
- use_case_implementation.R

#### R function (will be automatically loaded in the R script above):
- Fed_simulation_functions.R

#### Sample dataset:
- sample_data_to_run.csv

#### Instruction for replication of data application with sample dataset:
1. Download the R scripts under use case folder
2. Use R or RStudio to run the R script "use_case_implementation.R"
3. To generate the ROC curve as Figure 4 in the manuscript, please download and run the following corresponding R script:

   3.1. Figure4.R
