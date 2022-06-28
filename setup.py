#!/usr/bin/env python

from setuptools import setup, find_packages
import versioneer

NAME = "RIssmed"
DESCRIPTION = "The RNA Interaction via secondary structure mediation (RIssmed) tool suite analyses the change of RNA secondary structure upon binding other molecules."

# Set __version__ done by versioneer
# exec(open("NextSnakes/__init__.py").read())

setup(
    name=NAME,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Joerg Fallmann",
    author_email="fall@bioinf.uni-leipzig.de",
    scripts=[
        "RIssmed/ConstraintPLFold.py",
        "RIssmed/ConstraintFold.py",
        "RIssmed/CollectConsResults.py",
        "RIssmed/CollectWindowResults.py",
        "RIssmed/GenerateBigWig.py",
    ],
    packages=["RIssmed.Tweaks", "RIssmed", "RIssmed.vis"],
    package_data={"RIssmed.vis": ["assets/*", "templates/*"]},
    entry_points={
        "console_scripts": [
            "RIssmed_fold = RIssmed.ConstraintFold:main",
            "RIssmed_plfold = RIssmed.ConstraintPLFold:main",
            "RIssmed_collect_fold = RIssmed.CollectWindowResults:main",
            "RIssmed_collect_plfold = RIssmed.CollectConsResults:main",
            "RIssmed_visualize = RIssmed.visualize:main",
            "RIssmed_generate_bw = RIssmed.GenerateBigWig:main",
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
    url="https://github.com/jfallmann/RIssmed",
    description=DESCRIPTION,
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    zip_safe=False,
)
