import re
import pyopenms as oms

SCAN_PATTERN = r"(?:spectrum|scan)=(\d+)"


def extract_scan_id(spectrum: oms.MSSpectrum) -> str:
    """
    Extracts the scan ID from a given spectrum's native ID.

    Parameters
    ----------
    spectrum : oms.MSSpectrum
      The spectrum from which to extract the scan ID.

    Returns
    -------
    str
       The extracted scan ID if found, otherwise the original native ID.
    """
    match = re.search(SCAN_PATTERN, spectrum.getNativeID())
    if match:
        return match.group(1)
    return spectrum.getNativeID()
