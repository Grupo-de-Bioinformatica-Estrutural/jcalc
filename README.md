[![PyPI version](https://badge.fury.io/py/jcalc.svg)](https://badge.fury.io/py/jcalc)
[![Build Status](https://travis-ci.com/Joaodemeirelles/jcalc.svg?branch=main)](https://travis-ci.com/Joaodemeirelles/jcalc)
# jcalc
Python module to calculate vicinal coupling constant from Molecular Dynamics. A Jupyter Notebook with a [quickstart example](https://github.com/Joaodemeirelles/jcalc/blob/main/examples/Quickstart.ipynb) and [adding hydrogens example](https://github.com/Joaodemeirelles/jcalc/blob/main/examples/Adding_hydrogen.ipynb) shows how to use JCalc.

## Getting started

### Installing jcalc
* [Install GROMACS](http://www.gromacs.org/)
* [Install obabel](https://github.com/openbabel/openbabel)

### Install locally
```bash
pip install jcalc
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
