from distutils.core import setup
from setuptools import find_packages

package_name = 'm6anormalization'

setup(name=package_name,
      description="Calculate normalization constants for m6a levels",
      version='1.0',
      author="Alejandro A. Edera",
      author_email="edera.alejandro@gmail.com",
      packages=find_packages(),
      entry_points = {
          'console_scripts': [
              '{0} = {0}:main'.format(package_name)
          ]
      },
      setup_requires=['numpy'],
      install_requires=['numpy']
)
