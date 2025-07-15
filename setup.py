# setup.py
from setuptools import setup, find_packages

setup(
  name="flair_test_suite",
  version="0.1.0",
  packages=find_packages(where="src"),
  package_dir={"":"src"},
  install_requires=["toml"],
  entry_points={
    "console_scripts": [
      "flair-test-suite = flair_test_suite.cli:main"
    ]
  }
)
