# pyquantms
[![Python application](https://github.com/bigbio/pyquantms/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/pyquantms/actions/workflows/python-app.yml)
[![Python package](https://github.com/bigbio/pyquantms/actions/workflows/python-package.yml/badge.svg)](https://github.com/bigbio/pyquantms/actions/workflows/python-package.yml)
[![PyPI version](https://badge.fury.io/py/pyquantms.svg)](https://badge.fury.io/py/pyquantms)
[![Documentation Status](https://readthedocs.org/projects/pyquantms/badge/?version=latest)](https://pyquantms.readthedocs.io/en/latest/?badge=latest)

Python package with scripts and functions for the [quantms workflow](https://github.com/bigbio/quantms) for the analysis of quantitative proteomics data.

The package is available on PyPI: [pyquantms](https://pypi.org/project/pyquantms/)
```
pip install pyquantms
```

The following functionalities are available in the package:

### Diann scripts

- `dianncfg` - Create a configuration file for Diann including enzymes, modifications, and other parameters.
- `diann2mztab` - Convert Diann output to mzTab format. In addition, convert DIA-NN output to MSstats, Triqler or mzTab.
    The output formats are used for quality control and downstream analysis in quantms.

### SDRF scripts

- `openms2sample` - Extra sample information from OpenMS experimental design file. An example of OpenMS experimental design file is available [here](https://github.com/bigbio/pyquantms/blob/dev/tests/test_data/BSA_design_urls.tsv).
- `checksamplesheet` - Check the sample sheet for errors and inconsistencies. The experimental design coult be an OpenMS experimental design file or and SDRF file. 


