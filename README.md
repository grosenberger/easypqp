EasyPQP: Simple library generation for OpenSWATH
================================================

EasyPQP is a Python package that provides simplified and fast peptide query parameter generation for OpenSWATH. It can process input from MSFragger or other database search engines in pepXML format. Statistical validation can be conducted either using PyProphet or PeptideProphet/iProphet. Retention times are calibrated using an internal or external standard. In addition to a cumulative library, run-specific libraries are generated for non-linear RT alignment in OpenSWATH. For generation of PTM specific libraries that utilizes a unimod.xml database, you can further restrict the unimod.xml database file for modifications and site-specificities of interest.

Installation
============

We strongly advice to install EasyPQP in a Python [*virtualenv*](https://virtualenv.pypa.io/en/stable/). EasyPQP is compatible with Python 3.

Install the development version of *easypqp* from GitHub:

````
    $ pip install git+https://github.com/grosenberger/easypqp.git@master
````

Running EasyPQP
===============

*EasyPQP* is not only a Python package, but also a command line tool:

````
   $ easypqp --help
````

or:

````
   $ easypqp convert --help
   $ easypqp library --help
   $ easypqp reduce --help
   $ easypqp filter-unimod --help
````

Docker
======

EasyPQP is also available from Docker (automated builds):

Pull the development version of *easypqp* from DockerHub (synced with GitHub):

````
    $ docker pull grosenberger/easypqp:latest
````
