from setuptools import setup, find_packages

setup(
    name="LoReMiNE",
    version="1.0",
    description="Long Read-based Microbial genome mining pipeline",
    license='MIT',
    author="Amay Ajaykumar Agrawal",
    maintainer="Amay Ajaykumar Agrawal",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.11",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    package_dir={"src": "src"},           # <── tells setuptools code is under src/
    packages=["src"],
    install_requires=[
        "pandas==2.3.2",
        "biopython==1.81",
        "sortedcontainers==2.4.0",
        "sqlalchemy==2.0.2",
        "click==8.1.7",
    ],
    python_requires=">=3.11",
    keywords="bioinformatics",
    entry_points={
        "console_scripts": [
            "loremine=src.main:main",  # CLI entry point
        ],
    }
)