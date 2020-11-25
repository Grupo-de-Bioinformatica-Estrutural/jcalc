# *Running JCalc*


## Example (from pip)
### Calculate J values through Molecular Dynamics
    JCalcMD

### Calculate J values in structure (PDB format)
    JCalcPDB



## Example (from Docker)
### Calculate J values through Molecular Dynamics
    docker run -v $(pwd):/home/data jcalc:latest --x idose_sim.xtc \
    --t idose_sim.tpr --r HAID --suf 1 --skip 1000 --j j_input.txt \
    --ff add_hydrogen.ff

### Calculate J values in structure (PDB format)
    docker run -v $(pwd):/home/data jcalc:latest --p struct_file.pdb \
    --j j_input.txt --r HAID
