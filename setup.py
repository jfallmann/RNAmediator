#!/usr/bin/env python

from setuptools import setup, find_packages
import versioneer

NAME = "RNAmediator"
DESCRIPTION = "The RNA Interaction via secondary structure mediation (RNAmediator) tool suite analyses the change of RNA secondary structure upon binding other molecules."

# Set __version__ done by versioneer
# exec(open("NextSnakes/__init__.py").read())

setup(
    name=NAME,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Joerg Fallmann",
    author_email="fall@bioinf.uni-leipzig.de",
    scripts=[
        "RNAmediator/ConstraintPLFold.py",
        "RNAmediator/ConstraintFold.py",
        "RNAmediator/CollectConsResults.py",
        "RNAmediator/CollectWindowResults.py",
        "RNAmediator/GenerateBigWig.py",
    ],
    packages=["RNAmediator.Tweaks", "RNAmediator", "RNAmediator.vis"],
    package_data={"RNAmediator.vis": ["assets/*", "templates/*"]},
    entry_points={
        "console_scripts": [
            "RNAmediator_fold = RNAmediator.ConstraintFold:main",
            "RNAmediator_plfold = RNAmediator.ConstraintPLFold:main",
            "RNAmediator_collect_fold = RNAmediator.CollectWindowResults:main",
            "RNAmediator_collect_plfold = RNAmediator.CollectConsResults:main",
            "RNAmediator_visualize = RNAmediator.visualize:main",
            "RNAmediator_generate_bw = RNAmediator.GenerateBigWig:main",
        ]
    },
    include_package_data=True,
    install_requires=[
        "numpy",
        "biopython",
        "natsort",
        "pandas",
        "pybigwig",
        "dash",
        "dash-bootstrap-components",
    ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    license="LICENSE",
    url="https://github.com/jfallmann/RNAmediator",
    description=DESCRIPTION,
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    zip_safe=False,
)
