# CC2p-Event-Selection
The purpose of this code is to perform different portions of the event selection for the CC2p analysis. All of this code can be found in my app area on the GPVMS (/uboone/app/users/sfehlber/Analysis/CC2p). This code was developed using
uboonecode v08_00_00_52. To prepare to run anything from this repository, run the following lines:

```
  source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
  source /grid/fermiapp/products/uboone/setup_uboone.sh
  setup uboonecode v08_00_00_52 -q e17:prof
```

## PeLEE_ntuples Folder
This folder contains code to perform the CC2p event selectiion over the PeLEE ntuples for different MicroBooNE data products.

### Running the Event Selection
Inside of this folder you will find several scripts with the naming structure "twoproton_pelee" and then a sample identifier. Each pair of files runs the CC2p selection for that specific sample identifier's PELEE ntuples. 
You will not that therre is not much difference between the BNB, EXT, and dirt samples other than the application of a POT weight factor for the EXT and a POT weight and a MC weight applied for the dirt. 
The overlay sample is a bit more omplicated as an MC analysis is also performed. When running any of these pieces of code, a user prompt will appear asking which Run you wish to process. 
Type one of the corresponding numbers and then click enter. The script will run and produce an output root file (saved at /uboone/data/users/sfehlber/CC2p/PeLEE) and will print out various pieces of information to the screen such as the number
of events that remain after each selection cut, 

To run one of the scripts, do the following (after having setup you favorite version of uboonecode):
```
  root -b twoproton_pelee_IDENTIFIER.C
  twoproton_pelee_IDENTIFIER s
  s.Loop()
  USER INPUT PROMPT
```  
  In addition to these scripts there is a set of mc_efficiency scripts. These calculate the efficiency curves from the Overlay sample. It will also ask the user which run they want to process. To run this script, do the following
  (after having setup your favorite version of uboonecode):
  ```
    root -b mc_efficiency.C
    mc_efficiency s
    s.Loop()
    USER INPUT PROMPT
  ```
  ### Additional files
  
  - constants.h: Contains definitions the momentum thresholds, the track score cuts, the PID cut, various counters, and other constants. 
  - histograms_funcs.h: Class of functions for defining the histograms, filling them, and writing them.
  - variables.h: Class that contains functions to calculate different variables such as the detector angles, opening angles, the STVs, etc.
  - cuts.h/helper_funcs.h: Class of functions used to perform different cuts such as the FV check, check on the number of PFPs, and the reconstructed momentum.
  - neutrino_flux.py: Scrript to get the neutrino flux estimate. Run using python.
 
## NuWro Folder
This folder has a similar structure to the PeLEE_ntuples folder, but it looks at the NuWro samples that Afro prroduced. It is run in the exact same way as the above scripts, but with the indeitifier "nuwro". 

## Systematics Folder
Contains code to perform the CC2p selection over all the detector variation samples. Runs in exact same way as mentioned above. 

## Analysis Folder
The scripts in this folder are no longer up to do date. Please ignore this directory and instead consider the CC2p-Plotting-Tools repository.
