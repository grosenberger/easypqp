import pandas as pd
import pyopenms as po
import easypqp.openswathassaygenerator as OSAG


def main():

    infile = "/Users/owentsai/easypqp/test_library.tsv"
    in_type = None

    outfile = "/Users/owentsai/easypqp/easypqp/test_output.tsv"
    out_type = None ### if (in_type != None) else self._get_file_type(outfile)

    min_transitions = 6
    max_transitions = 6

    allowed_fragment_type = "b,y"
    allowed_fragment_charges = "1,2,3,4"

    enable_detection_specific_losses = False
    enable_detection_unspecific_losses = False

    precursor_mz_threshold = 0.025
    precursor_lower_mz_limit = 400.0
    precursor_upper_mz_limit = 1200.0
    product_mz_threshold = 0.025
    product_lower_mz_limit = 350.0
    product_upper_mz_limit = 2000.0

    swath_windows_file = None
    unimod_file = "/Users/owentsai/easypqp/easypqp/data/unimod.xml"
    enable_ipf = False
    max_num_alternative_localizations = 10000
    disable_identification_ms2_precursors = False
    disable_identification_specific_losses = False
    enable_identification_unspecific_losses = False
    enable_swath_specifity = False

    test_OSAG = OSAG.OpenSwathAssayGenerator(infile, in_type, outfile, out_type, min_transitions, max_transitions,allowed_fragment_type, allowed_fragment_charges, enable_detection_specific_losses,
    enable_detection_unspecific_losses, precursor_mz_threshold, precursor_lower_mz_limit, precursor_upper_mz_limit, product_mz_threshold, product_lower_mz_limit, product_upper_mz_limit, swath_windows_file, unimod_file, enable_ipf, max_num_alternative_localizations,
    disable_identification_ms2_precursors, disable_identification_specific_losses, enable_identification_unspecific_losses, enable_swath_specifity)

    
    test_OSAG.read_input_file()
    test_OSAG.annotate_transitions()
    test_OSAG.write_output_file()
    
    
if __name__ == '__main__':
    main()





