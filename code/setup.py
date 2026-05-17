from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="ubsumo-em-simulator",
    version="0.1.0",
    author="math-science agent",
    description="Ubiquitin/SUMO Electromagnetic Competition Simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
        "openmm>=8.0",
        "mdanalysis>=2.0",
        "h5py>=3.0",
    ],
    extras_require={
        "ml": ["deepmd-kit>=2.0"],
        "enhanced": ["plumed>=2.8"],
        "dev": ["pytest", "pytest-cov", "black", "mypy"],
    },
)