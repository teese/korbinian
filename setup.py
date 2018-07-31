from setuptools import setup, find_packages
from codecs import open
from os import path

# Get the long description from the readme file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'readme.rst')) as f:
    long_description = f.read()

classifiers = """\
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Medical Science Apps.
Topic :: Scientific/Engineering :: Chemistry
"""
		
setup(name="korbinian",
	author="Mark Teese",
	author_email="mark.teese@checkmytumhomepage.de",
	url="https://github.com/teese/korbinian",
	download_url = 'https://github.com/teese/korbinian/archive/0.3.4.tar.gz',
	description = "Bioinformatic sequence analysis of membrane proteins.",
	long_description=long_description,
	long_description_content_type='text/x-rst',
	packages=find_packages(),
	install_requires=['numpy','matplotlib','pandas', 'seaborn', 'scipy', 'biopython', 'tqdm'],
    license='MIT',
	project_urls=
    {
        'Wiki': 'https://github.com/teese/korbinian/wiki',
        'LangoschLab':'http://cbp.wzw.tum.de/index.php?id=9',
        "TU_Muenchen":"https://www.tum.de",
        "Freising":"https://en.wikipedia.org/wiki/Freising",
        "KorbinianSaint": "https://en.wikipedia.org/wiki/Corbinian",
        "KorbinianBeer" : "https://www.weihenstephaner.de/unsere-biere/korbinian"
    },
    include_package_data=True,
	keywords="bioinformatics protein sequence membrane conservation evolution membranous soluble polypeptide BLAST SIMAP",
	version="0.3.4"
	)
