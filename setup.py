from setuptools import find_packages, setup

VERSION = "0.0.3"

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
    "sdrf-pipelines==0.0.28",
    "pyopenms",
    "ms2rescore==3.0.2",
    "psm-utils==0.8.0",
    "pydantic",
    "pandas",
    "numpy",
    "pyarrow",
]
PYTHON_REQUIRES = ">=3.8,<4"

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name=NAME,
    version=VERSION,
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
