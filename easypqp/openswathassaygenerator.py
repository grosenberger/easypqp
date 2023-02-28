import os
import pandas as pd
import pyopenms as po
import ctypes
import click
from typing import Union, Tuple

from .targetedfileconverter import TargetedExperiment

def check_argument_values(arg_name: str, arg_value: any, expected_type: Tuple[Union[type, None], Union[Tuple, None]]) -> None:
    """
    Check if the given argument value is of the expected type and value range (if applicable).
    Raise a TypeError or ValueError if the value is invalid.
    """
    expected_type, expected_range = expected_type
    if isinstance(expected_type, list) and None in expected_type:
        pass
    elif not isinstance(arg_value, expected_type):
        raise TypeError(f"{arg_name} should be of type {expected_type.__name__} not type {arg_value.__class__}.")
    if expected_range is not None:
        # Handle numeric range
        if isinstance(expected_range, tuple) and len(expected_range) == 2:
            if not (expected_range[0] <= arg_value <= expected_range[1]):
                raise ValueError(f"{arg_name} should be within the range {expected_range}, cannot except {arg_value}.")
        elif isinstance(expected_range, list) and arg_value not in expected_range:
            raise ValueError(f"{arg_name} should be one of {expected_range}, cannot except '{arg_value}'.")

class OpenSwathAssayGenerator(TargetedExperiment):
    def __init__(self, 
                 infile: str, 
                 in_type: Union[str, None]=None, 
                 outfile: str="library.pqp", 
                 out_type: Union[str, None]=None,
                 min_transitions: int=6,
                 max_transitions: int=6,
                 allowed_fragment_type: str="b,y",
                 allowed_fragment_charges: str='1,2,3,4',
                 enable_detection_specific_losses: bool=False,
                 enable_detection_unspecific_losses: bool=False,
                 precursor_mz_threshold: float=0.025,
                 precursor_lower_mz_limit: float=400.0,
                 precursor_upper_mz_limit: float=1200.0,
                 product_mz_threshold: float=0.025,
                 product_lower_mz_limit: float=350.0,
                 product_upper_mz_limit: float=2000.0,
                 swath_windows_file: str="",
                 unimod_file: str="",
                 enable_ipf: bool=False,
                 max_num_alternative_localizations: int=10000,
                 disable_identification_ms2_precursors: bool=False,
                 disable_identification_specific_losses: bool=False,
                 enable_identification_unspecific_losses: bool=False,
                 enable_swath_specifity: bool=False) -> None:
        super().__init__(True)

        # Valdiate arguments
        check_argument_values("infile", infile, (str, None))
        check_argument_values("outfile", outfile, (str, None))
        # Handle types
        if in_type is None:
            in_type = self._get_file_type(infile)
        if out_type is None:
            out_type = self._get_file_type(outfile)
        check_argument_values("in_type", in_type, ([str, None], ['tsv', 'mrm', 'pqp', 'TraML']))
        check_argument_values("out_type", out_type, ([str, None], ['tsv', 'pqp', 'TraML']))
        check_argument_values("product_mz_threshold", product_mz_threshold, (float, None))
        check_argument_values("allowed_fragment_types", allowed_fragment_types, (str, None)) # TODO: Add value check, to ensure valid fragment types
        check_argument_values("allowed_fragment_charges", allowed_fragment_charges, (str, None)) # TODO: Add value check to ensure ints are in string of charges
        check_argument_values("enable_detection_specific_losses", enable_detection_specific_losses, (bool, None))
        check_argument_values("enable_detection_unspecific_losses", enable_detection_unspecific_losses, (bool, None))


    def generate_openswathassays():
        return None