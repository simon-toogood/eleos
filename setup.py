"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name="nemesis_eleos",  # Required
    version="0.1",  # Required
    description="A Python interface to the NEMESIS spectral inversion tool",  # Optional
    long_description=long_description,  # Optional
    long_description_content_type="text/markdown",  # Optional ()
    url="https://github.com/simon-toogood/eleos",  # Optional
    author="Simon Toogood",  # Optional
    author_email="scat2@leicester.ac.uk",  # Optional
    classifiers=[  # Optional
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="astronomy, atmospheric physics, spectral inversion, NEMESIS",  # Optional
    packages=["nemesis_eleos"],  # Required
    python_requires=">=3.10, <4",
    # If there are data files included in your packages that need to be
    # installed, specify them here.
    include_package_data=True,
    package_dir={"": "."},
    # package_data={  # Optional
    #     "nemesis_eleos": ["data/*"],
    # },
    # Entry points. The following would provide a command called `sample` which
    # executes the function `main` from this package when invoked:
    #entry_points={  # Optional
    #    "console_scripts": [
    #        "sample=sample:main",
    #    ],
    #},
    project_urls={  # Optional
        "Source": "https://github.com/simon-toogood/eleos/",
    },
)