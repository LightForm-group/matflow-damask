"""Pip installation script."""

import os
import re
from setuptools import find_packages, setup


def get_version():

    ver_file = 'matflow_damask/_version.py'
    with open(ver_file) as handle:
        ver_str_line = handle.read()

    ver_pattern = r'^__version__ = [\'"]([^\'"]*)[\'"]'
    match = re.search(ver_pattern, ver_str_line, re.M)
    if match:
        ver_str = match.group(1)
    else:
        msg = 'Unable to find version string in "{}"'.format(ver_file)
        raise RuntimeError(msg)

    return ver_str


def get_long_description():
    readme_file = 'README.md'
    with open(readme_file, encoding='utf-8') as handle:
        contents = handle.read()
    return contents


def get_changelog():
    changelog_file = 'CHANGELOG.md'
    with open(changelog_file, encoding='utf-8') as handle:
        contents = handle.read()
    return contents


setup(
    author="Adam J. Plowman, Michael D. Atkinson, Guy L. P. Bowker",
    author_email='adam.plowman@manchester.ac.uk',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Matflow extension for DAMASK.",
    entry_points="""
        [matflow.extension]
        damask=matflow_damask
    """,
    install_requires=[
        'matflow',
        'vtk',
        'damask',
        'damask-parse',
    ],
    license="MIT license",
    long_description=get_long_description() + '\n\n' + get_changelog(),
    long_description_content_type='text/markdown',
    keywords='matflow, materials-science, computational-workflow',
    name='matflow-damask',
    packages=find_packages(),
    project_urls={
        'GitHub': 'https://github.com/LightForm-group/matflow-damask'
    },
    version=get_version(),
)
