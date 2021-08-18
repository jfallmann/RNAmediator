from setuptools import setup, find_packages


setup(
    name="RIssmed",
    version="0.0.1",
    author="Joerg Fallmann",
    author_email="fall@bioinf.uni-leipzig.de",
    scripts=[
        "RIssmed/ConstraintPLFold.py",
        "RIssmed/ConstraintFold.py",
        "RIssmed/CollectConsResults.py",
        "RIssmed/CollectWindowResults.py",
    ],
    packages=["RIssmed.Tweaks", "RIssmed"],
    entry_points={
        "console_scripts": [
            "RIssmed_fold = RIssmed.ConstraintFold:main",
            "RIssmed_plfold = RIssmed.ConstraintPLFold:main",
            "RIssmed_collect_fold = RIssmed.CollectWindowResults:main",
            "RIssmed_collect_plfold = RIssmed.CollectConsResults:main",
            "RIssmed_visualize = RIssmed.visualize:main",
        ]
    },
    license="LICENSE.txt",
    url="https://github.com/jfallmann/RIssmed",
    long_description_content_type="text/markdown",
    description="The RNA Interaction via secondary structure mediation (RIssmed) tool suite analyses the change of RNA "
    "secondary structure upon binding other molecules. (Most of the functionality currently missing)",
    long_description=open(
        "RIssmed/Tweaks/README.md"
    ).read(),  # Change this to include other RISSmed stuff
    include_package_data=True,
    install_requires=["numpy"],
)
