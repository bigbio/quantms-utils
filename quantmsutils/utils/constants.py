"""
 These columns are now aligned with the column names in quantms.io definition. Here is the quantms.io definition:
 https://github.com/bigbio/quantms.io
"""

CHARGE = "charge"  # charge of the precursor ion
SCAN = "scan"  # Definition in quantms.io https://github.com/bigbio/quantms.io/blob/main/docs/README.adoc#scan
MS_LEVEL = "ms_level"
NUM_PEAKS = "num_peaks"  # number of peaks in the spectrum (MS1 or MS2)
BASE_PEAK_INTENSITY = "base_peak_intensity"
SUMMED_PEAK_INTENSITY = "summed_peak_intensities"
MAX_INTENSITY = "max_intensity"
RETENTION_TIME = "rt"
EXPERIMENTAL_MASS_TO_CHARGE = "observed_mz"
ACQUISITION_DATETIME = "acquisition_datetime"
MONOISOTOPIC_MZ = "monoisotopic_mz"

MZ_ARRAY = "mz_array"
INTENSITY_ARRAY = "intensity_array"
