from setuptools import setup, find_packages
import os

# Parse version string
this_directory = os.path.dirname(os.path.abspath(__file__))
version_file = os.path.join(this_directory, "psix", "_version_.py")
exec(open(version_file).read())


setup(
    name="psix",
    version=__version__,
    packages=find_packages(),

    install_requires=[
        'matplotlib>=3.2.1',
        'numba>=0.43.1',
        'numpy>=1.18.1',
        'seaborn>=0.10.0',
        'scipy>=1.4.1',
        'pandas>=1.0.3',
        'tqdm>=4.46.0',
        'statsmodels>=0.11.0',
        'scikit-learn>=0.22.2',
        'anndata>=0.7.5',
    ],

    include_package_data=True,

    author="Carlos F Buen Abad Najar",
    author_email="cfbuenabadn@berkeley.edu",
    description="Psix is a computational tool that finds cell-state associated exons in scRNA-seq data",
    keywords="",
    url="",
    license=""
)
