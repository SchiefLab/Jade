#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0',
    # TODO: Put package requirements here
]

setup_requirements = [
    # TODO(jadolfbr): Put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    # TODO: Put package test requirements here
]

setup(
    name='jade',
    version='1',
    description="A repository for modules and applications to aid in the design and analysis of Biological molecules, especially when working with Rosetta or PyRosetta",
    long_description=readme + '\n\n' + history,
    author="Jared Adolf-Bryfogle",
    author_email='jadolfbr@gmail.com',
    url='https://github.com/jadolfbr/jade',
    packages=find_packages(include=['jade']),
    entry_points={
        'console_scripts': [
            'jade=jade.cli:main',
        ],
    },
    include_package_data=True,
    install_requires=requirements,
    license="BSD license",
    zip_safe=False,
    keywords='jade',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
