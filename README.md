# jcalc
Python module to calculate vicinal coupling constant from Molecular Dynamics

## Getting started

### Installing jcalc
* [Install GROMACS](http://www.gromacs.org/)
* [Install obabel](https://github.com/openbabel/openbabel)

### Install locally
```bash
./easy_install.py
```

### Install by docker
```bash
docker pull jlmeirelles/jcalc:latest
```

## Running jcalc

### Running on simulation (XTC and TPR file)
```bash
jcalc -x sim.xtc -t sim.tpr -n j_input.tsv
```

### Running on directory with steps as PDB files
```bash
jcalc -i pdb_dir -n j_input.tsv
```

### Running on single PDB file
```bash
jcalc -p file.pdb -i j_input.tsv
```

### Running on docker
```bash
docker run -v $(pwd):/home/data jlmeirelles/jcalc -x sim.xtc \
-t sim.tpr -n j_input.tsv
```
