#!/usr/bin/env python3

from setuptools import setup

with open('README.rst') as f:
    README = f.read()

setup(name='pyespresso',
      version='0.0.1',
      classifiers=[
          'Development Status :: 0.0.1',
          'License :: MIT License',
          'Programming Language :: Python :: 3.5',
      ],
      description='The Pythonic Quantum ESPRESSO package',
      long_description=README,
      url='https://bitbucket.org/singularitti/pyespresso',
      author='Qi Zhang',
      author_email='qz2280@columbia.edu',
      license='MIT',
      install_requires=['addict',
                        'beeprint',
                        'json_tricks',
                        'lazy-property',
                        'matplotlib',
                        'numpy',
                        'scipy',
                        'spglib'
                        ],
      include_package_data=True,
      packages=['basics',
                'calculators',
                'miscellaneous',
                'plotters',
                'readers',
                'submitters'
                ],
      zip_safe=False)
