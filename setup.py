import sys
from pathlib import Path

from setuptools import find_packages, setup

try:
    from pasmopy import __author__, __email__, __maintainer__
except ImportError:
    __author__ = __maintainer__ = "Hiroaki Imoto"
    __email__ = "himoto@protein.osaka-u.ac.jp"


def get_version() -> str:
    """Read version from file"""
    version_filepath = Path(__file__).parent / "pasmopy" / "version.py"
    with open(version_filepath) as f:
        for line in f:
            if line.startswith("__version__"):
                return line.strip().split()[-1][1:-1]
    assert False


def setup_package():
    # Python version check.
    if sys.version_info[:2] < (3, 7):
        raise RuntimeError("Pasmopy requires at least Python version 3.7")

    setup(
        name="pasmopy",
        version=get_version(),
        description="Patient-Specific Modeling in Python",
        long_description=Path("README.md").read_text("utf-8"),
        long_description_content_type="text/markdown",
        license="Apache 2.0",
        author=__author__,
        author_email=__email__,
        maintainer=__maintainer__,
        maintainer_email=__email__,
        url="https://github.com/pasmopy/pasmopy",
        download_url="https://github.com/pasmopy/pasmopy/releases",
        project_urls={
            "Documentation": "https://pasmopy.readthedocs.io/en/latest/",
            "Source Code": "https://github.com/pasmopy/pasmopy",
            "Bug Tracker": "https://github.com/pasmopy/pasmopy/issues",
        },
        packages=find_packages(exclude=["tests", "docs"]),
        install_requires=[
            l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
        ],
        extras_require={
            "dev": [
                "black==22.1.0",
                "flake8==4.0.1",
                "isort==5.10.1",
                "pre-commit",
                "pytest",
            ],
            "docs": [
                "sphinx>=1.7",
                "sphinx_rtd_theme>=0.3",
                "sphinx_autodoc_typehints>=1.10",
                "sphinxcontrib-bibtex>=2.2",
            ],
        },
        python_requires=">=3.7",
        keywords=[
            "systems",
            "biology",
            "cancer",
            "stratification",
            "personalized",
            "modeling",
            "simulation",
            "precision",
            "oncology",
        ],
        classifiers=[
            "Intended Audience :: Science/Research",
            "Intended Audience :: Healthcare Industry",
            "License :: OSI Approved :: Apache Software License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Medical Science Apps.",
            "Topic :: Software Development",
            "Topic :: Software Development :: Libraries",
            "Topic :: Software Development :: Libraries :: Python Modules",
        ],
        zip_safe=False,
    )


if __name__ == "__main__":
    setup_package()
