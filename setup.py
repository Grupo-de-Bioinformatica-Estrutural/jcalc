# Base imports
from setuptools import setup

# Project imports
from jcalc.settings import JCALC_VERSION

setup(
    name="jcalc",
    version=JCALC_VERSION,
    description="jcalc: Calculate NMR J values from \
        Molecular Dynamics simulations",
    url="https://github.com/Joaodemeirelles/jcalc/",
    author="Joao Luiz de Meirelles",
    author_email="jldemeirelles@gmail.com",
    licence="Academic",
    packages=[
        "jcalc",
        "jcalc.core",
        "jcalc.logger",
    ],
    scripts=["jcalc/jcalc"],
    install_requires=[
        "numpy==1.19.4",
        "statistics==1.0.3.5",
        "biopython==1.78",
        "argparse==1.4.0",
        "pathlib==1.0.1",
        "flake8==3.9.0"
    ],
    zip_safe=False
)
