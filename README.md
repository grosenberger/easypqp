EasyPQP: Simple library generation for OpenSWATH
================================================

[![CI](https://github.com/grosenberger/easypqp/actions/workflows/ci.yml/badge.svg)](https://github.com/grosenberger/easypqp/actions/workflows/ci.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/easypqp?link=https%3A%2F%2Fpypi.org%2Fproject%2Feasypqp%2F)](https://pypi.org/project/easypqp/)

EasyPQP is a Python package that provides simplified and fast peptide query parameter generation for OpenSWATH. It can process input from MSFragger, Sage or other database search engines in pepXML/idXML/tsv format. Statistical validation can be conducted either using PyProphet or PeptideProphet/iProphet. Retention times and ion mobilities are calibrated using an internal or external standard. In addition to a cumulative library, run-specific libraries are generated for non-linear RT alignment in OpenSWATH. For generation of PTM specific libraries that utilizes a unimod.xml database, you can further restrict the unimod.xml database file for modifications and site-specificities of interest. It also supports in-silico library generation.

Installation
============

We strongly advice to install EasyPQP in a Python [*virtualenv*](https://virtualenv.pypa.io/en/stable/). EasyPQP is compatible with Python 3.

Install the development version of *easypqp* from GitHub:

````
    $ pip install git+https://github.com/grosenberger/easypqp.git@master
````

### Full Installation

To install all optional features:

````
    $ pip install easypqp[all]
````

This will install the `easypqp_rs` package, which provides the in-silico library generation feature and pyprophet for statistical validation.

Running EasyPQP
===============

*EasyPQP* is not only a Python package, but also a command line tool:

````
   $ easypqp --help
````

or:

````
   $ easypqp convert --help
   $ easypqp convertpsm --help
   $ easypqp convertsage --help
   $ easypqp library --help
   $ easypqp insilico-library --help
   $ easypqp reduce --help
   $ easypqp filter-unimod --help
   $ easypqp openswath-assay-generator --help
   $ easypqp openswath-decoy-generator --help
   $ easypqp targeted-file-converter --help
````

Generating an *In-Silico* Library
=================================

The in-silico library generation feature is included if you installed EasyPQP with the `[all]`  or `[rust]` extras (to install the `easypqp_rs` package).

To generate an in-silico library, you can use the `insilico-library` command. For example:

````
   $ easypqp insilico-library --fasta your_proteome.fasta --output_file insilico_library.tsv
````

For more information on the parameters and JSON configuration file, see the [Configuration Reference](https://github.com/singjc/easypqp-rs?tab=readme-ov-file#configuration-reference)

> [!NOTE]
> If no `retention_time`, `ion_mobility`, or `ms2_intensity` fields are provided under `dl_feature_generators` in the config, pretrained models will be automatically downloaded and used. The current default pretrained models used are:
> - RT: `rt_cnn_tf` - A CNN-Transformer model trained on the [ProteomicsML repository RT dataset](https://proteomicsml.org/datasets/retentiontime/ProteomeTools_RT.html). This model is based on AlphaPeptDeep's CNN-LSTM implementation, with the biLSTM replaced by a Transformer encoder.
> - CCS: `ccs_cnn_tf` - A CNN-Transformer model trained on the [ProteomicsML repository CCS dataset](https://proteomicsml.org/datasets/ionmobility/Meier_TIMS.html). This model is also based on AlphaPeptDeep's CNN-LSTM implementation, with the biLSTM replaced by a Transformer encoder.
> - MS2: `ms2_bert` - A BERT-based model retreived from AlphaPeptDeep's pretrained models.

If you want just a standalone portable rust binary, you can download one from the [easypqp-rs releases page](https://github.com/singjc/easypqp-rs/releases).

Docker
======

EasyPQP is also available from Docker (automated builds):

Pull the development version of *easypqp* from DockerHub (synced with GitHub):

````
    $ docker pull grosenberger/easypqp:latest
````
