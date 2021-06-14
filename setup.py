from setuptools import setup, find_packages

from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='easypqp',
    version='0.1.14',
    description='EasyPQP: Simple library generation for OpenSWATH',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/grosenberger/easypqp",
	author="George Rosenberger",
	author_email="gr2578@cumc.columbia.edu",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    include_package_data=True,
    install_requires=['Click','numpy','scipy','scikit-learn','statsmodels','pandas>=1.1.0','biopython','pyopenms>=2.6.0','matplotlib','seaborn'],
    # extras_require={  # Optional
    #     'dev': ['check-manifest'],
    #     'test': ['coverage'],
    # },
    package_data={  # Optional
        'easypqp': ['data/unimod.xml'],
    },
    entry_points={
        'console_scripts': [
            'easypqp=easypqp.main:cli',
        ],
    },
    extras_require={
        'PyProphet': ['pyprophet'],
    },
)
