from setuptools import setup, find_packages

setup(name="korbinian_sw_mp",
      author="Mark Teese",
      packages=find_packages(),
      install_requires=['numpy','matplotlib','pandas', 'biopython'],
      version="0.2.9")