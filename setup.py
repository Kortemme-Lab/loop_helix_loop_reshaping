#!/usr/bin/env python2

from distutils.core import setup

setup(
    name='loop_helix_loop_reshaping',
    version='0.0.0',
    author='Xingjie Pan',
    author_email='xingjiepan@gmail.com',
    url='https://github.com/xingjiepan/loop_helix_loop_reshaping',
    packages=[
        'loop_helix_loop_reshaping',
    ],
    install_requires=[
        'numpy',
        'matplotlib',
        'docopt',
    ],
    extras_require = {
        'weblogo':  ['weblogo'],
    },
    entry_points={
        'console_scripts': [
        ],
    },
    description='PyRosetta scripts for reshaping a part of a protein into a loop-helix-loop unit',
    long_description=open('README.md').read(),
    classifiers=[
        'Programming Language :: Python :: 2',
        'Intended Audience :: Science/Research',
    ],
)
