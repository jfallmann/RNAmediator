from setuptools import setup, find_packages
import versioneer

NAME = "RIssmed"
DESCRIPTION = "The RNA Interaction via secondary structure mediation (RIssmed) tool suite analyses the change of RNA secondary structure upon binding other molecules. (Most of the functionality currently missing)"

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
    license="LICENSE",
    url="https://github.com/jfallmann/RIssmed",
    description=DESCRIPTION,
    long_description=open(
        "RIssmed/Tweaks/README.md"
    ).read(),  # Change this to include other RISSmed stuff
    long_description_content_type="text/markdown",
    include_package_data=True,
    install_requires=["numpy"],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
)
