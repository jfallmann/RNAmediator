from setuptools import setup, find_packages
import versioneer

setup(
    name="RIssmed",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
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
            "RIssmed fold = RIssmed.ConstraintFold:main",
            "RIssmed plfold = RIssmed.ConstraintPLFold:main",
            "RIssmed collect fold = RIssmed.CollectWindowResults:main",
            "RIssmed collect plfold = RIssmed.CollectConsResults:main",
            "RIssmed visualize = RIssmed.visualize:main",
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
