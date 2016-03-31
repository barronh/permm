# PERMM Installation #

Known Platforms: 
 - Unix
 - Mac
 - Linux

Requirements:
 - [Python >=2.5](http://python.org)
 - [NumPy >=1.2](http://numpy.scipy.org)
 - [matplotlib =0.98.3](http://matplotlib.sourceforge.net)
 - [YAML](http://www.yaml.org)
 - [PseudoNetCDF >=1](https://dawes.sph.unc.edu/trac/PseudoNetCDF)
 - One of the below netcdf readers
  - [pynetcdf](http://pypi.python.org/pypi/pynetcdf/0.7)
  - [netCDF3](http://code.google.com/p/netcdf4-python)
  - [pupynere](http://pypi.python.org/pypi/pupynere/)
  - [Scientific](http://dirac.cnrs-orleans.fr/plone/software/scientificpython/)

Instructions
 1. [Download current source distribution](https://github.com/barronh/permm/archive/master.zip)
 1. unzip permm-master.zip -d permm-master
 1. cd permm-master/trunk
 1. python setup.py install