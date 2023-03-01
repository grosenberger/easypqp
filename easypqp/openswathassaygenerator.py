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

def check_fragment_type(input_str: str):
    possible_fragment_types = ['b','y','a','x','c','z']
    if input_str not in possible_fragment_types:
        raise ValueError(f"{input_str} is not one of the possible fragment types {possible_fragment_types}")

def string_to_list(input_str: str, output_type: type):
    ### TODO
    return None


def read_swath_file():
    ### TODO
    return None


class OpenSwathAssayGenerator(TargetedExperiment):
    def __init__(self, infile, in_type, outfile, out_type, min_transitions, max_transitions, allowed_fragment_type, allowed_fragment_charges, enable_detection_specific_losses, enable_detection_unspecific_losses, precursor_mz_threshold, precursor_lower_mz_limit,
                 precursor_upper_mz_limit, product_mz_threshold, product_lower_mz_limit, product_upper_mz_limit, swath_windows_file, unimod_file, enable_ipf, max_num_alternative_localizations, disable_identification_ms2_precursors, disable_identification_specific_losses, enable_identification_unspecific_losses, enable_swath_specifity) -> None:
        super().__init__(True)

        self.infile = infile 
        self.in_type = in_type

        self.outfile = outfile
        self.out_type = out_type ### if (in_type != None) else self._get_file_type(outfile)

        self.min_transitions = min_transitions
        self.max_transitions = max_transitions

        allowed_fragment_type_string = allowed_fragment_type.split(",")
        allowed_fragment_types_bytes = []
        for s in allowed_fragment_type_string:
            allowed_fragment_types_bytes.append(bytes(s, encoding='utf-8')) #utf-8
        self.allowed_fragment_type = allowed_fragment_types_bytes
        print(self.allowed_fragment_type)

        allowed_fragment_charges_string = allowed_fragment_charges.split(",")
        allowed_fragment_types_int = []
        for s in allowed_fragment_charges_string:
            allowed_fragment_types_int.append(int(s))
        self.allowed_fragment_charges = allowed_fragment_types_int
        print(self.allowed_fragment_charges)

        self.enable_detection_specific_losses = enable_detection_specific_losses
        self.enable_detection_unspecific_losses = enable_detection_unspecific_losses
        self.precursor_mz_threshold = precursor_mz_threshold
        self.precursor_lower_mz_limit = precursor_lower_mz_limit
        self.precursor_upper_mz_limit = precursor_upper_mz_limit
        self.product_mz_threshold = product_mz_threshold
        self.product_lower_mz_limit = product_lower_mz_limit
        self.product_upper_mz_limit = product_upper_mz_limit

        ### TODO: read input swath file
        self.swath_windows_file = swath_windows_file
        ### TODO: read unimod file 
        self.unimod_file = unimod_file
        ### TODO: implement enable ipf
        self.enable_ipf = enable_ipf
        self.max_num_alternative_localizations = max_num_alternative_localizations
        self.disable_identification_ms2_precursors = disable_identification_ms2_precursors
        self.disable_identification_specific_losses = disable_identification_specific_losses
        self.enable_identification_unspecific_losses = enable_identification_unspecific_losses
        self.enable_swath_specifity = enable_swath_specifity



        ### check argument
        # # Valdiate arguments
        # check_argument_values("infile", infile, (str, None))
        # check_argument_values("outfile", outfile, (str, None))
        # # Handle types
        # if in_type is None:
        #     in_type = self._get_file_type(infile)
        # if out_type is None:
        #     out_type = self._get_file_type(outfile)
        # check_argument_values("in_type", in_type, ([str, None], ['tsv', 'mrm', 'pqp', 'TraML']))
        # check_argument_values("out_type", out_type, ([str, None], ['tsv', 'pqp', 'TraML']))
        # check_argument_values("product_mz_threshold", product_mz_threshold, (float, None))
        # check_argument_values("allowed_fragment_types", allowed_fragment_types, (str, None)) # TODO: Add value check, to ensure valid fragment types
        # check_argument_values("allowed_fragment_charges", allowed_fragment_charges, (str, None)) # TODO: Add value check to ensure ints are in string of charges
        # check_argument_values("enable_detection_specific_losses", enable_detection_specific_losses, (bool, None))
        # check_argument_values("enable_detection_unspecific_losses", enable_detection_unspecific_losses, (bool, None))

    def read_input_file(self) -> None:
        self.load_library(self.infile, self.in_type)

    def annotate_transitions(self) -> None:
        click.echo("Info: Annotating transitions")
        assays = po.MRMAssay()
        assays.reannotateTransitions(self.tr_exp, self.precursor_mz_threshold, self.product_mz_threshold, self.allowed_fragment_type, self.allowed_fragment_charges, self.enable_detection_specific_losses, self.enable_detection_unspecific_losses, -4) ### todo convert fragment type to bytes

        click.echo("Info: Annotating detecting transitions")
        assays.restrictTransitions(self.tr_exp, self.product_lower_mz_limit, self.product_upper_mz_limit, list(list())) #TODO read swath file and save 
        assays.detectingTransitions(self.tr_exp, self.min_transitions, self.max_transitions)

    def write_output_file(self) -> None:
        self.write_library(self.outfile, self.out_type)
