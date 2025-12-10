from datetime import datetime
import json
import tempfile
from typing import Union
import click


def timestamped_echo(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    click.echo(f"{timestamp} - {message}")




def create_json_config(as_bytes: bool = False) -> Union[str, bytes]:
    """
    Create a JSON configuration file for EasyPQP In-silico library generation.
    """
    config = {
        "version": "0.1.0",
        "database": {
            "enzyme": {
                "missed_cleavages": 1,
                "min_len": None,
                "max_len": None,
                "cleave_at": "KR",
                "restrict": "P",
                "c_terminal": None,
                "semi_enzymatic": None
            },
            "peptide_min_mass": 500.0,
            "peptide_max_mass": 5000.0,
            "static_mods": {
                "C": 57.0215
            },
            "variable_mods": {},
            "max_variable_mods": 2,
            "decoy_tag": "rev_",
            "generate_decoys": True,
            "fasta": ""
        },
        "insilico_settings": {
            "precursor_charge": [2, 4],
            "max_fragment_charge": 1,
            "min_transitions": 6,
            "max_transitions": 6,
            "fragmentation_model": "cid",
            "allowed_fragment_types": ["b", "y"],
            "rt_scale": 100.0
        },
        "dl_feature_generators": {
            "device": "cpu",
            "fine_tune_config": {
                "fine_tune": False,
                "train_data_path": "",
                "batch_size": 256,
                "epochs": 3,
                "learning_rate": 0.001,
                "save_model": True
            },
            "instrument": "QE",
            "nce": 20.0,
            "batch_size": 64
        },
        "peptide_chunking": 0,
        "output_file": "./easypqp_insilico_library.tsv",
        "write_report": True,
        "parquet_output": False
    }

    json_str = json.dumps(config, indent=2)

    if as_bytes:
        return json_str.encode('utf-8')
    else:
        with tempfile.NamedTemporaryFile('w+', suffix=".json", delete=False) as tmp:
            tmp.write(json_str)
            tmp.flush()
            return tmp.name
