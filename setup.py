from setuptools import find_packages, setup
import os
import codecs


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")


NAME = "quantms-utils"
LICENSE = "MIT License"
DESCRIPTION = "Python package with scripts and helpers for the QuantMS workflow"
AUTHOR = "Yasset Perez-Riverol, Dai Chengxin"
AUTHOR_EMAIL = "ypriverol@gmail.com"
URL = "https://www.github.com/bigbio/pyquantms"
PROJECT_URLS = {
    "Documentation": "https://docs.quantms.org/en/latest/",
    "quantms Workflow": "https://github.com/bigbio/quantms",
    "Tracker": "https://github.com/bigbio/pyquantms/issues",
}

KEYWORDS = [
    "quantms",
    "Proteomics",
]

CLASSIFIERS = [
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Apache Software License",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python :: 3 :: Only",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Development Status :: 5 - Production/Stable",
]

INSTALL_REQUIRES = [
    "click",
    "sdrf-pipelines>=0.0.29",
    "pyopenms",
    "ms2rescore==3.0.3",
    "deeplc==2.2.38",
    "ms2pip==4.0.0.dev8",
    "psm-utils==0.8.2",
    "deeplcretrainer",
    "pydantic",
    "pandas",
    "protobuf>=3.9.2,< 4",
    "numpy",
    "pyarrow",
    "pygam",
    "scipy==1.13.1",
    "scikit-learn",
]
PYTHON_REQUIRES = ">=3.8,<4"

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name=NAME,
    version=get_version("quantmsutils/__init__.py"),
    license=LICENSE,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    project_urls=PROJECT_URLS,
    keywords=KEYWORDS,
    classifiers=CLASSIFIERS,
    packages=find_packages(),
    include_package_data=True,
    entry_points={"console_scripts": ["quantmsutilsc=quantmsutils.quantmsutilsc:main"]},
    install_requires=INSTALL_REQUIRES,
    python_requires=PYTHON_REQUIRES,
)
