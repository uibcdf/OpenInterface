from setuptools import setup, find_packages
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
#import distutils.extension

extensions_list=[]

setup(
    name='openinterface',
    version='0.0.0',
    author='UIBCDF Lab',
    author_email='uibcdf@gmail.com',
    package_dir={'openinterface': 'openinterface'},
    packages=find_packages(),
    ext_modules=extensions_list,
    package_data={'openinterface': []},
    scripts=[],
    url='http://uibcdf.org',
    download_url ='https://github.com/uibcdf/OpenInterface',
    license='MIT',
    description="---",
    long_description="---",
)
