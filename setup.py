from setuptools import setup, find_packages

setup(
    name="msmsrt_ssvm_dataset",
    version="0.0.1",
    license="MIT",
    packages=find_packages(exclude=["metfrag", "tests", "examples", "*.ipynb"]),

    # Minimum requirements the package was tested with
    install_requires=[
        "numpy",
        "scipy",
        "scikit-learn",
        "pandas",
        "matplotlib",
        "seaborn",
        "joblib",
        "xlrd",
        "setuptools>=46.1"
    ],

    # Metadata
    author="Eric Bach",
    author_email="eric.bach@aalto.fi",
    description="TODO",
    url="https://github.com/bachi55/massbank_for_ML",
)