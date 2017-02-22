from setuptools import setup, find_packages

setup(name="korbinian",
      author="Mark Teese",
      packages=find_packages(),
      install_requires=['numpy','matplotlib','pandas', 'biopython'],
      version="0.3.0")