import pyopenms as po
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

class OpenSwathDecoyGenerator(TargetedExperiment):
    def __init__(self, 
                 infile: str, 
                 outfile: str="library.pqp", 
                 in_type: Union[str, None]=None, 
                 out_type: Union[str, None]=None, 
                 method: str="shuffle",
                 decoy_tag: str="DECOY_",
                 min_decoy_fraction: float=0.8,
                 aim_decoy_fraction: float=1.0,
                 shuffle_max_attempts: int=30,
                 shuffle_sequence_identity_threshold: float=0.5,
                 shift_precursor_mz_shift: float=0.0,
                 shift_product_mz_shift: float=20.0,
                 product_mz_threshold: float=0.025,
                 allowed_fragment_types: str="b,y",
                 allowed_fragment_charges: str="1,2,3,4",
                 enable_detection_specific_losses: bool=False,
                 enable_detection_unspecific_losses: bool=False,
                 switchKR: bool=True,
                 separate: bool=False) -> None:
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
        check_argument_values("method", method, (str, ['shuffle', 'pseudo-reverse', 'reverse', 'shift']))
        check_argument_values("decoy_tag", decoy_tag, (str, None))
        check_argument_values("min_decoy_fraction", min_decoy_fraction, (float, (0, 1)))
        check_argument_values("aim_decoy_fraction", aim_decoy_fraction, (float, (0, 1)))
        check_argument_values("shuffle_max_attempts", shuffle_max_attempts, (int, None))
        check_argument_values("shuffle_sequence_identity_threshold", shuffle_sequence_identity_threshold, (float, (0, 1)))
        check_argument_values("shift_precursor_mz_shift", shift_precursor_mz_shift, (float, None))
        check_argument_values("shift_product_mz_shift", shift_product_mz_shift, (float, None))
        check_argument_values("product_mz_threshold", product_mz_threshold, (float, None))
        check_argument_values("allowed_fragment_types", allowed_fragment_types, (str, None)) # TODO: Add value check, to ensure valid fragment types
        check_argument_values("allowed_fragment_charges", allowed_fragment_charges, (str, None)) # TODO: Add value check to ensure ints are in string of charges
        check_argument_values("enable_detection_specific_losses", enable_detection_specific_losses, (bool, None))
        check_argument_values("enable_detection_unspecific_losses", enable_detection_unspecific_losses, (bool, None))
        check_argument_values("switchKR", switchKR, (bool, None))
        check_argument_values("separate", separate, (bool, None))

        # TODO: Move this up before argument validation for specific arg?
        # Transform string
        allowed_fragment_types = allowed_fragment_types.split(",")
        allowed_fragment_types = [s.encode('utf-8') for s in allowed_fragment_types]
        allowed_fragment_charges = allowed_fragment_charges.split(",")
        allowed_fragment_charges = [int(charge) for charge in allowed_fragment_charges]

        # Assign values to self
        for name, value in locals().items():
            if name != 'self':
                # print(f"Info: Setting {name} = {value}")
                setattr(self, name, value)

        # Load target experiment
        self.load_library(self.infile, self.in_type)

    def generate_decoys(self) -> None:
        # Initiate decoy experiment
        self.tr_decoy = po.TargetedExperiment()

        # Generate decoys
        decoys = po.MRMDecoy()
        decoys.generateDecoys(self.tr_exp, self.tr_decoy, self.method, self.aim_decoy_fraction, self.switchKR, self.decoy_tag, self.shuffle_max_attempts, self.shuffle_sequence_identity_threshold, self.shift_precursor_mz_shift, self.shift_product_mz_shift, self.product_mz_threshold, self.allowed_fragment_types, self.allowed_fragment_charges, self.enable_detection_specific_losses, self.enable_detection_unspecific_losses, -4)

        click.echo(f"Info: Number of target peptides: {len(self.tr_exp.getPeptides())}")
        click.echo(f"Info: Number of decoy peptides: {len(self.tr_decoy.getPeptides())}")
        click.echo(f"Info: Number of target proteins: {len(self.tr_exp.getProteins())}")
        click.echo(f"Info: Number of decoy proteins: {len(self.tr_decoy.getProteins())}")
        
        if len(self.tr_decoy.getPeptides()) / len(self.tr_exp.getPeptides()) < self.min_decoy_fraction or len(self.tr_decoy.getProteins()) / len(self.tr_exp.getProteins()) < self.min_decoy_fraction:
            raise click.ClickException(f"The number of decoys for peptides or proteins is below the threshold of {(self.min_decoy_fraction * 100)}% of the number of targets.")
        
        if self.separate:
            click.echo(f"Info: Writing only decoys to file: {self.outfile}")
            self.tr_exp = self.tr_decoy
        else:
            click.echo(f"Info: Writing targets and decoys to file: {self.outfile}")
            self.tr_exp += self.tr_decoy

        self.write_library(self.outfile, self.out_type)