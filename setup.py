from setuptools import find_packages, setup

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="RNAHyperFold",
    packages=find_packages(where="."),
    version="0.1.0",
    description="A library that provides a number of classes and methods designed to facilitate the analysis and visualization of changes in RNA folding at different temperatures using hypergraph-based approaches.",
    author="Nicol Buratti",
    install_requires=requirements,
    setup_requires=[],
    tests_require=[],
    test_suite="",
)
