# quantms-utils
[![Python application](https://github.com/bigbio/quantms-utils/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/quantms-utils/actions/workflows/python-app.yml)
[![Python package](https://github.com/bigbio/quantms-utils/actions/workflows/python-package.yml/badge.svg)](https://github.com/bigbio/quantms-utils/actions/workflows/python-package.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/ea6903630b3a4d15b674a16b8ce594a7)](https://app.codacy.com/gh/bigbio/quantms-utils/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![PyPI version](https://badge.fury.io/py/quantms-utils.svg)](https://badge.fury.io/py/quantms-utils)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Python package with scripts and functions for the [quantms workflow](https://github.com/bigbio/quantms) for the analysis of quantitative proteomics data.

The package is available on PyPI: [quantms-utils](https://pypi.org/project/quantms-utils/)
```
pip install quantms-utils
```

## Available Scripts 

The following functionalities are available in the package:

### Diann scripts

- `dianncfg` - Create a configuration file for Diann including enzymes, modifications, and other parameters.
- `diann2mztab` - Convert Diann output to mzTab format. In addition, convert DIA-NN output to MSstats, Triqler or mzTab.
    The output formats are used for quality control and downstream analysis in quantms.

### SDRF scripts

- `openms2sample` - Extra sample information from OpenMS experimental design file. An example of OpenMS experimental design file is available [here](https://github.com/bigbio/quantms-utils/blob/dev/tests/test_data/BSA_design_urls.tsv).
- `checksamplesheet` - Check the sample sheet for errors and inconsistencies. The experimental design coult be an OpenMS experimental design file or and SDRF file. 

### Other scripts

- `psmconvert` - The convert_psm function converts peptide spectrum matches (PSMs) from an idXML file to a CSV file, optionally filtering out decoy matches. It extracts and processes data from both the idXML and an associated spectra file, handling multiple search engines and scoring systems.
- `mzmlstats` - The `mzmlstats` processes mass spectrometry data files in either `.mzML` or `Bruker .d` formats to extract and compile statistics about the spectra. It supports generating detailed or ID-only CSV files based on the spectra data.

#### mzml statistics 

quantms-utils have multiple scripts to generate mzML stats. These files are used by multiple tools and packages within quantms ecosystem for quality control, mzTab generation, etc. Here are some details about the formats, the fields they contain and gow they are computed.

<details>
<summary>MS info and details</summary>

`mzmlstats` allows the user to produce a file containing all features for every signal in the MS/MS experiment. The produced file is a parquet file, with the original name of the file plus the following postfix `{file_name}_ms_info.parquet`. Here, the definition of each column and how they are estimated and used: 

- `scan`: The scan accession for each MS and MS/MS signal in the mzML, depending on the manufacturer, the scan will have different formats. Example, for thermo (e.g `controllerType=0 controllerNumber=1 scan=43920`). We tried to find the definition of [quantms.io](https://github.com/bigbio/quantms.io/blob/main/docs/README.adoc#scan). 
- `ms_level`: The MS level of the signal, 1 for MS and 2 for MS/MS.
- `num_peaks`: The number of peaks in the MS. Compute with pyopenms with `spectrum.get_peaks()`.
- `base_peak_intensity`: The max intensity in the spectrum (MS or MS/MS).
- `summed_peak_intensities`: The sum of all intensities in the spectrum (MS or MS/MS).
- `rt`: The retention time of the spectrum, capture with pyopenms with `spectrum.getRT()`.

For MS/MS signals, we have the following additional columns:

- `precursor_charge`: The charge of the precursor ion, if the signal is MS/MS. Capture with pyopenms with `spectrum.getPrecursors()[0].getCharge()`.
- `precursor_mz`: The m/z of the precursor ion, if the signal is MS/MS. Capture with pyopenms with `spectrum.getPrecursors()[0].getMZ()`.
- `precursor_intensity`: The intensity of the precursor ion, if the signal is MS/MS. Capture with pyopenms with `spectrum.getPrecursors()[0].getIntensity()`. If the precursor is not annotated (present), we use the purity object to get the information; see note below. 
- `precursor_rt`: The retention time of the precursor ion, if the signal is MS/MS. See note below.
- `precursor_total_intensity`: The total intensity of the precursor ion, if the signal is MS/MS. See note below.
 
> **NOTE**: For all the precursor-related information, we are using the first precursor in the spectrum. The following columns `intensity` (if not annotated), `precursor_rt`, and `precursor_total_intensity` we use the following pyopnems code: 
> ```python
> precursor_spectrum = mzml_exp.getSpectrum(precursor_spectrum_index)
> precursor_rt = precursor_spectrum.getRT()
> purity = oms.PrecursorPurity().computePrecursorPurity(precursor_spectrum, precursor, 100, True)
> precursor_intensity = purity.target_intensity
> total_intensity = purity.total_intensity
> ```

</details>

<details>
<summary>MS2 info and details</summary>

`mzmlstats` allows the user to produce a file containing all the MS2 spectra including the intesities and masses of every peak. The produced file is a parquet file, with the original name of the file plus the following postfix `{file_name}_ms2_info.parquet`. Here, the definition of each column and how they are estimated and used:

- `scan`: The scan accession for each MS and MS/MS signal in the mzML, depending on the manufacturer, the scan will have different formats. Example, for thermo (e.g `controllerType=0 controllerNumber=1 scan=43920`). We tried to find the definition of [quantms.io](https://github.com/bigbio/quantms.io/blob/main/docs/README.adoc#scan).
- `ms_level`: The MS level of the signal, all of them will be 2.
- `mz_array`: The m/z array of the peaks in the MS/MS signal. Capture with pyopenms with `mz_array, intensity_array = spectrum.get_peaks()`.
- `intensity_array`: The intensity array of the peaks in the MS/MS signal. Capture with pyopenms with `mz_array, intensity_array = spectrum.get_peaks()`.

</details>

<details>
<summary>MS1 Feature Maps</summary>

We use the `FeatureFinderMultiplexAlgorithm` from [OpenMS](https://pyopenms.readthedocs.io/en/latest/apidocs/_autosummary/pyopenms/pyopenms.FeatureFinderMultiplexAlgorithm.html) 
to extract the features from the MS1 spectra. We use an algorithm based on the original implementation by [Andy Lin](https://doi.org/10.1093/bioinformatics/btad058). The output of this algorithm is a feature map, which contains the following information:

- `feature_mz`: The m/z of the feature.
- `feature_rt`: The retention time of the feature.
- `feature_intensity`: The intensity of the feature.
- `feature_charge`: The charge of the feature.
- `feature_quality`: The quality of the feature.
- `feature_percentile_tic`: The percentile of the feature in the total ion current.
- `feature_id`: The unique identifier of the feature generated by OpenMS. 
- `feature_min_rt`: The minimum retention time of the feature within the feature map.
- `feature_min_mz`: The minimum m/z of the feature within the feature map.
- `feature_max_rt`: The maximum retention time of the feature within the feature map.
- `feature_max_mz`: The maximum m/z of the feature within the feature map.
- `feature_num_scans`: The number of scans that the feature is present in the feature map.
- `feature_scans`: The scans where the feature is present in the feature map.

The tool will generate a gzip compressed parquet file with the extension `{file_name}_ms1_feature_info.parquet`.

## Contributions and issues

Contributions and issues are welcome. Please, open an issue in the [GitHub repository](https://github.com/bigbio/quantms) or PR in the [GitHub repository](https://github.com/bigbio/quantms-utils).
