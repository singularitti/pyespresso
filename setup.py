#!/usr/bin/env python3

from markdown import markdownFromFile
from setuptools import setup

setup(name='PyQuE',
      version='0.0.1',
      classifiers=[
          'Development Status :: 0.0.1',
          'License :: MIT License',
          'Programming Language :: Python :: 3.5',
      ],
      description='The Pythonic Quantum ESPRESSO package',
      long_description=markdownFromFile(input='README.md'),
      url='https://bitbucket.org/singularitti/pyque',
      author='Qi Zhang',
      author_email='qz2280@columbia.edu',
      license='MIT',
      install_requires=['addict',
                        'beeprint',
                        'fortranformat',
                        'json_tricks',
                        'lazy-property',
                        'Markdown',
                        'matplotlib',
                        'numpy',
                        'scipy',
                        'spglib'],
      include_package_data=True,
      packages=['basics',
                'calculators',
                'miscellaneous',
                'plotters',
                'readers',
                'submitters'],
      zip_safe=False)
