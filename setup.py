from setuptools import find_packages
from numpy.distutils.core import Extension
import os

src_dir = 'snapping/Fortran'
ext1 = Extension('snapping',
                 [os.path.join(src_dir, 'addtau.f'),
                  os.path.join(src_dir, 'rollsort.f'),
                  os.path.join(src_dir, 'smerge.f'),
                  os.path.join(src_dir, 'sbeam.f'),
                  os.path.join(src_dir, 'rnext.f'),
                  os.path.join(src_dir, 'ascnx.f'),
                  os.path.join(src_dir, 'clump.f'),
                  os.path.join(src_dir, 'asc3d.f'),
                  os.path.join(src_dir, 'snapy.f')],
                 )

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'snapping',
          description       = "Impulse data beamforming",
          author            = "Edward Richards",
          author_email      = "edwardlrichards@ucsd.edu",
          version="0.1",
          packages=find_packages(),
          ext_modules = [ext1]
          )
