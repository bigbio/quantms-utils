"""
 MS1 Feature Detection, this algorithm is based on the OpenMS library and is used to detect MS1 features from mzML files.
 In addition, it adds some normalization and filtering steps from the
 previous algorithm by Andy Lim https://github.com/bmx8177/MS1Connect
 published in https://doi.org/10.1093/bioinformatics/btad058.

 This algorithm is used to detect MS1 features from mzML files and save them to parquet format.
"""

import bisect
import logging
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any, Union
import pandas as pd
from pyopenms import MzMLFile, MSExperiment, FeatureMap, LogType, PeakFileOptions
import pyopenms as oms

logging.basicConfig(format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


class MS1FeatureDetector:
    """
    Class for detecting MS1 features from mzML files and saving to parquet format.
    """

    def __init__(self, min_ptic: float = 0.05, max_ptic: float = 0.95, ms_level: int = 1):
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

        self.min_ptic = min_ptic
        self.max_ptic = max_ptic
        self.ms_level = ms_level

        # Initialize options for file loading
        self.options = PeakFileOptions()
        self.options.setMSLevels([self.ms_level])

    def _calc_tic(self, experiment: MSExperiment) -> float:
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

    def _get_ptic_data(self, experiment: MSExperiment) -> Tuple[List[float], List[float]]:
        """
        Convert TIC to pTIC (percentile TIC) for all MS scans. The pTIC is the cumulative sum of TIC
        up to a given retention time, divided by the total TIC.

        Parameters
        ----------
        experiment : MSExperiment
            The MS experiment containing spectra.

        Returns
        -------
        Tuple[List[float], List[float]]
            Lists of retention times and corresponding pTIC values.
        """
        total_tic = self._calc_tic(experiment)
        if total_tic == 0:
            logger.warning("Total TIC is zero, check input data")
            return [], []

        rt_list = []
        ptic_list = []
        sum_tic = 0.0

        logger.info("Converting TIC to pTIC")
        for scan in experiment:
            if scan.getMSLevel() == self.ms_level:
                mz, intensities = scan.get_peaks()
                rt_list.append(scan.getRT())
                ptic_list.append(sum_tic / total_tic)
                sum_tic += sum(intensities)

        return rt_list, ptic_list

    def _find_ptic_for_rt(self, rt: float, rt_list: List[float], ptic_list: List[float]) -> float:
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
        self, features: FeatureMap, rt_list: List[float], ptic_list: List[float]
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
            help(feature)
            mz = round(feature.getMZ(), 4)
            intensity = feature.getIntensity()
            rt = round(feature.getRT(), 4)
            charge = feature.getCharge()

            # Get quality metrics
            quality = feature.getOverallQuality()

            # Interpolate pTIC for this retention time
            ptic = round(self._find_ptic_for_rt(rt, rt_list, ptic_list), 4)

            # Filter by pTIC
            if self.min_ptic <= ptic <= self.max_ptic:
                feature_list.append(
                    {
                        "mz": mz,
                        "intensity": intensity,
                        "rt": rt,
                        "pTIC": ptic,
                        "charge": charge,
                        "quality": quality,
                        "id": feature.getUniqueId(),
                    }
                )
            else:
                logger.debug(f"Skipping feature at RT {rt} due to pTIC {ptic}")

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
            Path to output parquet file.
        top_n : int, optional
            Top N most intense features to save, by default 1000
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
            file_handler = MzMLFile()
            file_handler.setOptions(self.options)

            input_map = MSExperiment()
            file_handler.load(str(input_path), input_map)

            if input_map.size() == 0:
                logger.error("No spectra found in input file")
                return None

            # Get pTIC data
            rt_list, ptic_list = self._get_ptic_data(input_map)

            # Prepare for feature finding
            input_map.updateRanges()

            # Run feature finder
            logger.info("Running feature detection")
            feature_finder = oms.FeatureFinderAlgorithmPicked()
            features = oms.FeatureMap()
            seeds = oms.FeatureMap()
            params = (
                feature_finder.getParameters()
            )  # In the original implementation it was feature_finder.getParameters("centroided")
            feature_finder.run(input_map, features, params, seeds)
            features.setUniqueIds()

            logger.info(f"Found {features.size()} features")

            # Extract features
            feature_list = self._extract_features(features, rt_list, ptic_list)

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


# Example usage
if __name__ == "__main__":
    import argparse

    # Setup argument parser
    parser = argparse.ArgumentParser(description="MS1 Feature Detection")
    parser.add_argument("--input", "-i", required=True, help="Input mzML file or directory")
    parser.add_argument("--output", "-o", required=True, help="Output parquet file or directory")
    parser.add_argument("--top", "-t", type=int, default=1000, help="Top N features")
    parser.add_argument("--min-ptic", type=float, default=0.05, help="Minimum pTIC")
    parser.add_argument("--max-ptic", type=float, default=0.95, help="Maximum pTIC")

    args = parser.parse_args()
    # Create detector
    detector = MS1FeatureDetector(min_ptic=args.min_ptic, max_ptic=args.max_ptic)

    result = detector.process_file(args.input, args.output, args.top)
    if result:
        print(f"Successfully processed {args.input} to {result}")
    else:
        print(f"Failed to process {args.input}")
