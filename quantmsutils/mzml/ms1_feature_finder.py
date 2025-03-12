"""
 MS1 Feature Detection, this algorithm is based on the OpenMS library and is used to detect MS1 features from mzML files.
 In addition, it adds some normalization and filtering steps from the
 previous algorithm by Andy Lim https://github.com/bmx8177/MS1Connect
 published in https://doi.org/10.1093/bioinformatics/btad058.

 We improved the original algorithm with the following ideas:
  - Using FeatureFinderMultiplexAlgorithm instead of FeatureFinder as originally implemented by Andy Lim. This will provide
    a more robust way to perform FeatureFinding.
  - Remove the filtering of percentile TIC for features, we leave this step to future consuming tools of the data to perform
    extra curation of the features based on percentile_tic, or quality of the feature, etc.
  - We annotated additional features such as min and max retention time and mz values.
  - This algorithm is used to detect MS1 features from mzML files and save them to parquet format.
"""

import bisect
import logging
from pathlib import Path
from typing import List, Optional, Dict, Any, Union
import pandas as pd
import pyopenms as oms

from quantmsutils.openms import extract_scan_id

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


class MS1FeatureDetector:
    """
    Class for detecting MS1 features from mzML files and saving to parquet format.
    """

    def __init__(self, ms_level: int = 1):
        """
        Initialize the MS1 feature detector.

        Parameters
        ----------
        min_ptic : float, optional
            Minimum percentile TIC to include features, by default 0.05
        max_ptic : float, optional
            Maximum percentile TIC to include features, by default 0.95
        ms_level : int, optional
            MS level to analyze, by default 1
        """
        # Configure logging

        self.ms_level = ms_level
        # Initialize options for file loading
        self.options = oms.PeakFileOptions()
        self.options.setMSLevels([self.ms_level])

    def _calc_tic(self, experiment: oms.MSExperiment) -> float:
        """
        Calculate all the TIC in the MS experiment, summing all intensities from all scans.

        Parameters
        ----------
        experiment : MSExperiment
            The MS experiment containing spectra.

        Returns
        -------
        float
            Total ion current.
        """
        return sum(
            sum(intensity)
            for mz, intensity in (
                scan.get_peaks() for scan in experiment if scan.getMSLevel() == self.ms_level
            )
        )

    def _get_ptic_data(self, experiment: oms.MSExperiment):
        """
        Convert TIC to pTIC (percentile TIC) for all MS scans. The pTIC is the cumulative sum of TIC
        up to a given retention time, divided by the total TIC.

        Parameters
        ----------
        experiment : MSExperiment
            The MS experiment containing spectra.
        """
        total_tic = self._calc_tic(experiment)
        if total_tic == 0:
            logger.warning("Total TIC is zero, check input data")
            return [], []

        rt_list, ptic_list, scans = [], [], []
        sum_tic = 0.0

        logger.info("Converting TIC to pTIC")
        for scan in experiment:
            if scan.getMSLevel() == self.ms_level:
                intensities = scan.get_peaks()[1]
                rt_list.append(scan.getRT())
                ptic_list.append(sum_tic / total_tic)
                sum_tic += sum(intensities)
                scans.append(extract_scan_id(scan))

        return rt_list, ptic_list, scans

    @staticmethod
    def _find_ptic_for_rt(rt: float, rt_list: List[float], ptic_list: List[float]) -> float:
        """
        Find the pTIC value for a given retention time by interpolation.

        Parameters
        ----------
        rt : float
            Retention time to find pTIC for.
        rt_list : List[float]
            List of retention times.
        ptic_list : List[float]
            List of pTIC values corresponding to rt_list.

        Returns
        -------
        float
            Interpolated pTIC value.
        """
        if not rt_list or not ptic_list:
            return 0.0

        index = bisect.bisect_left(rt_list, rt)

        # Handle edge cases
        if index == 0:
            return ptic_list[0]
        if index >= len(rt_list):
            return ptic_list[-1]

        # Linear interpolation between adjacent points
        rt_left = rt_list[index - 1]
        rt_right = rt_list[index]
        ptic_left = ptic_list[index - 1]
        ptic_right = ptic_list[index]

        # Calculate interpolated pTIC
        rt_frac = (rt - rt_left) / (rt_right - rt_left) if rt_right != rt_left else 0
        return ptic_left + rt_frac * (ptic_right - ptic_left)

    def _extract_features(
        self,
        features: oms.FeatureMap,
        rt_list: List[float],
        ptic_list: List[float],
        scans: List[str],
    ) -> List[Dict[str, Any]]:
        """
        Extract feature information and filter by pTIC.

        Parameters
        ----------
        features : FeatureMap
            Feature map from feature finder.
        rt_list : List[float]
            List of retention times.
        ptic_list : List[float]
            List of pTIC values corresponding to rt_list.

        Returns
        -------
        List[Dict[str, Any]]
            List of feature dictionaries.
        """
        feature_list = []

        for feature in features:
            mz = round(feature.getMZ(), 4)
            intensity = feature.getIntensity()
            rt = round(feature.getRT(), 4)
            charge = feature.getCharge()

            # Get quality metrics
            quality = feature.getOverallQuality()

            # Interpolate pTIC for this retention time
            ptic = round(self._find_ptic_for_rt(rt, rt_list, ptic_list), 4)

            convexhull = feature.getConvexHull().getBoundingBox()

            minRT, minMZ = convexhull.minPosition()
            maxRT, maxMZ = convexhull.maxPosition()
            select_scans = self._get_selected_scans(scans, rt_list, minRT, maxRT)
            num_scans = len(select_scans)

            feature_list.append(
                {
                    "feature_mz": mz,
                    "feature_intensity": intensity,
                    "feature_rt": rt,
                    "feature_charge": charge,
                    "feature_percentile_tic": ptic,
                    "feature_quality": quality,
                    "feature_id": feature.getUniqueId(),
                    "feature_min_rt": minRT,
                    "feature_min_mz": minMZ,
                    "feature_max_rt": maxRT,
                    "feature_max_mz": maxMZ,
                    "feature_num_scans": num_scans,
                    "feature_scans": select_scans,
                }
            )

        return feature_list

    def process_file(
        self,
        input_file: Union[str, Path],
        output_file: Union[str, Path],
        sort_by: str = "intensity",
        ascending: bool = False,
    ) -> Optional[str]:
        """
        Process an mzML file to detect features and save to parquet format.

        Parameters
        ----------
        input_file : str
            Path to mzML file to process.
        output_file : str
            Path to an output parquet file.
        sort_by : str, optional
            Column to sort by, by default "intensity"
        ascending : bool, optional
            Sort order, by default False (descending)

        Returns
        -------
        Optional[str]
            Path to the output file if successful, None otherwise.
        """
        try:
            # Validate inputs
            input_path = Path(input_file)
            if not input_path.exists():
                logger.error(f"Input file not found: {input_file}")
                return None

            # Create output directory if needed
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Load MS data
            logger.info(f"Loading data from {input_file}")
            file_handler = oms.MzMLFile()
            experiment = oms.MSExperiment()
            picker = oms.PeakPickerHiRes()

            filtered_experiment = oms.MSExperiment()

            file_handler.setOptions(self.options)
            file_handler.load(str(input_path), experiment)
            picker.pickExperiment(experiment, filtered_experiment, False)
            experiment = filtered_experiment

            filtered_experiment = oms.MSExperiment()
            for spec in experiment:
                if spec.getMSLevel() == self.ms_level:
                    filtered_experiment.addSpectrum(spec)
            experiment = filtered_experiment

            if experiment.size() == 0:
                logger.error("No spectra found in input file")
                return None

            # Get pTIC data
            rt_list, ptic_list, scans = self._get_ptic_data(experiment)

            # Run feature finder
            logger.info("Running feature detection")
            feature_finder = oms.FeatureFinderMultiplexAlgorithm()
            params = feature_finder.getParameters()
            params.setValue("algorithm:labels", "[]")
            feature_finder.setParameters(params)

            feature_finder.run(exp=experiment, progress=True)
            features = feature_finder.getFeatureMap()
            features.setUniqueIds()

            logger.info(f"Found {features.size()} features")

            # Extract features
            feature_list = self._extract_features(features, rt_list, ptic_list, scans)

            # Create DataFrame
            df = pd.DataFrame(feature_list)

            # Sort and limit features
            if sort_by in df.columns:
                df = df.sort_values(by=sort_by, ascending=ascending)

            # Save to Parquet
            df.to_parquet(output_path, index=False, compression="gzip")
            logger.info(f"Saved {len(df)} features to {output_file}")

            return str(output_path)

        except Exception as e:
            logger.exception(f"Error processing {input_file}: {str(e)}")
            return None

    @staticmethod
    def _get_selected_scans(scans: List[str], rt_list: List[float], min_rt: float, max_rt: float):
        """
        This function returns the scans that are within the RT range of the feature.
        The scans and rt_list are two lists that have the same length and the same order.
        :param scans: List of scans ids
        :param rt_list: List of retention times
        :param min_rt: Minimum retention time for the feature
        :param max_rt: Maximum retention time for the feature
        :return:
        """
        selected_scans = []
        for i, rt in enumerate(rt_list):
            if min_rt <= rt <= max_rt:
                selected_scans.append(scans[i])
        return selected_scans
