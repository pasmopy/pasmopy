import os
import sys

from setuptools import find_packages, setup


def get_version() -> str:
    """Read version from file"""
    version_filepath = os.path.join(os.path.dirname(__file__), "dyaus", "version.py")
    with open(version_filepath) as f:
        for line in f:
            if line.startswith("__version__"):
                return line.strip().split()[-1][1:-1]


def main():
    # Python version check.
    if sys.version_info[:2] < (3, 7):
        sys.exit("Dyaus requires at least Python version 3.7")

    # set long_description and requirements
    here = os.path.abspath(os.path.dirname(__file__))
    long_description = open(os.path.join(here, "README.md")).read()
    requirements = open(os.path.join(here, "requirements.txt")).read()

    setup(
        name="dyaus",
        version=get_version(),
        description="Dynamics-driven automatic subtyping",
        long_description=long_description,
        long_description_content_type="text/markdown",
        license="Apache 2.0",
        author="Hiroaki Imoto",
        author_email="himoto@protein.osaka-u.ac.jp",
        url="https://github.com/okadalabipr/dyaus",
        packages=find_packages(),
        install_requires=requirements.splitlines(),
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
            "modeling",
        ],
        classifiers=[
            "Intended Audience :: Science/Research",
            "Intended Audience :: Healthcare Industry",
            "License :: OSI Approved :: Apache Software License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3 :: Only",
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
