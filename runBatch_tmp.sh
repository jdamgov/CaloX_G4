#!/bin/bash
#$-V
#$-cwd
#$-S /bin/bash
#$-N g4cry_A0
#$-o logs/$JOB_NAME.stdout$JOB_ID
#$-e logs/$JOB_NAME.stderr$JOB_ID
#$-q omni
#$-P quanah

# module load gnu geant4 root

# set G4 and ROOT environment
source ~/geant4.10.05.p01-install/share/Geant4-10.5.1/geant4make/geant4make.sh
source /usr/local/root/bin/thisroot.sh
source ./muonSetupMac.sh

export CaloXG4OutName="XXXXXYYYYY/sampl_Si_0p5_XXXXXYYYYY_ZZZZZ.root"  # output root tree file name.
export CaloXG4RunNumber=1   # run number
export CaloXG4EventNumber=1 # initial event number
export CaloXG4SeedA=1       # 0=default (fixed), 1 auto, >1 manual seed setting
export CaloXG4SeedB=2019    # second seed for manual seed setting (CaloXG4SeedA>1)
export CaloXG4SaveSeeds=1   # 0=not save,  1=save seeds. 

export CaloXG4map1CellEdepON=1 # Hit Map 0=off, 1=on
export CaloXG4map2PidEdepON=1  # Edep per particle (Pid) 0=off, 1=on
export CaloXG4mapTcutMin=0       # time cut on Map creation (minimum)  (psec)  
export CaloXG4mapTcutMax=100000  # time cut on Map creation (maximum) (psec)  
export CaloXG4map1CelEcutMinSave=10 # energy cut on Edep per cell during tree output (kev).

# run a job
/home/jdamgov/CaloX/sim/exampleB4b -m XXXXXYYYYY/paramBatch.mac

