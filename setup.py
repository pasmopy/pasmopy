import sys
from pathlib import Path

from setuptools import find_packages, setup

try:
    from dyaus import __author__, __email__, __maintainer__, __version__
except ImportError:
    __author__ = ", ".join(["Hiroaki Imoto", "Sawa Yamashiro"])
    __maintainer__ = "Hiroaki Imoto"
    __email__ = "himoto@protein.osaka-u.ac.jp"
    __version__ = "0.0.1"


def main():
    # Python version check.
    if sys.version_info[:2] < (3, 7):
        sys.exit("Dyaus requires at least Python version 3.7")

    setup(
        name="dyaus",
        version=__version__,
        description="Dynamics-driven automatic subtyping",
        long_description=Path("README.md").read_text("utf-8"),
        long_description_content_type="text/markdown",
        license="Apache 2.0",
        author=__author__,
        author_email=__email__,
        maintainer=__maintainer__,
        maintainer_email=__email__,
        url="https://github.com/dyaus-dev/dyaus",
        packages=find_packages(exclude=["tests"]),
        install_requires=[
            l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
        ],
        extras_require={
            "dev": [
                "black==20.8b1",
                "flake8",
                "isort",
                "pre-commit",
                "pytest",
            ]
        },
        python_requires=">=3.7",
        keywords=[
            "systems",
            "biology",
            "cancer",
            "subtype",
            "classification",
            "modeling",
            "simulation",
            "transcriptome",
            "ccle",
            "tcga",
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
            "Programming Language :: R",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Medical Science Apps.",
            "Topic :: Software Development",
            "Topic :: Software Development :: Libraries",
            "Topic :: Software Development :: Libraries :: Python Modules",
        ],
    )


if __name__ == "__main__":
    main()
