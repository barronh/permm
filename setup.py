from distutils.core import setup
try:
    from setuptools import setup
except:
    pass
import os
import sys
from warnings import warn

def find_packages():
    import os
    packages = []
    walker = os.walk('src')
    prefix = os.path.join(os.path.curdir,'src')
    for thisdir, itsdirs, itsfiles in walker:
        if '__init__.py' in itsfiles:
            packages.append(thisdir[len(prefix)-1:])
    
    return packages
            
def find_data():
    import os
    import re
    data_pattern = re.compile(r'.*(.|_)(yaml|nc|net|irr|phy|ptb|sum|voc|txt|xls|graffle)$')
    data = []
    prefix = os.path.join(os.path.curdir,'src', 'permm')
    walker = os.walk('src')
    for thisdir, itsdirs, itsfiles in walker:
        if thisdir != os.path.join('src','permm.egg-info'):
            data.extend([os.path.join(thisdir[len(prefix)-1:],f) for f in itsfiles if data_pattern.match(f) is not None])
    
    return data

packages = find_packages()
data = find_data()

setup(name = 'permm',
      version = '1.0',
      author = 'Barron Henderson',
      author_email = 'barronh@gmail.com',
      maintainer = 'Barron Henderson',
      maintainer_email = 'barronh@gmail.com',
      description = 'Python Environment for Reaction Mechanisms/Mathematics provides dynamic analysis tools for evaluating chemical networks easily.',
      packages = packages,
      package_dir = {'': 'src'},
      package_data = {'permm': data},
      scripts = ['scripts/permm'],
      requires = ['numpy (>=1.2)', 'yaml', 'netCDF4'],
      url = 'http://github.com/barronh/permm/',
      download_url = 'https://github.com/barronh/permm/archive/v1.0.zip'
      )
