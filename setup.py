from distutils.core import setup
try:
    from setuptools import setup
except:
    pass
import os
import sys
from warnings import warn
netcdf_pkgs = [('netCDF4', 'Dataset', 'Variable'), \
               ('netCDF3', 'Dataset', 'Variable'), \
               ('pupynere', 'netcdf_file', 'netcdf_variable')]
for pkg, reader, var in netcdf_pkgs:
    try:
        NetCDFFile = getattr(__import__(pkg, fromlist = [reader]),reader)
        print >> file(os.path.join('src', 'permm', 'netcdf.py'),'wb'), """
__all__ = ['NetCDFFile']
__doc__ = \"\"\"
.. _netcdf
:mod:`netcdf` -- netcdf import point
====================================

.. module:: netcdf
   :platform: Unix, Windows
   :synopsis: Povides a single import point for a package.  If
              a user has one of many netcdf interfaces, this module
              selects it and provides it.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
\"\"\"
from %s import %s as NetCDFFile
from %s import %s as NetCDFVariable
""" % (pkg, reader, pkg, var)
        break
    except ImportError, e:
        warn(e.message)
else:
    raise ImportError, "Did not find a NetCDFFile object"

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
      version = '1.0rc',
      author = 'Barron Henderson',
      author_email = 'barronh@gmail.com',
      maintainer = 'Barron Henderson',
      maintainer_email = 'barronh@gmail.com',
      packages = packages,
      package_dir = {'': 'src'},
      package_data = {'permm': data},
      requires = [pkg, 'numpy (>=1.2)', 'yaml']
      )
