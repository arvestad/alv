import setuptools
import sys 

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('alv/version.py') as fh:
    exec(fh.read())
    
if sys.version_info.major < 3:
    sys.exit('\n'
             'Sorry, Python 2 is not supported\n'
             'Did you run pip install alv?\n'
             'Try \'pip3 install alv\'')

elif sys.version_info.minor < 2:
    sys.exit('\nSorry, Python < 3.2 is not supported\n')
    
setuptools.setup(
    name="alv",
    version=__version__,
    author="Lars Arvestad",
    author_email="arve@math.su.se",
    description="A console-based sequence alignment viewer",
    scripts = ['bin/alv'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/arvestad/alv",
    test_suite = "tests",
    packages=setuptools.find_packages(),
    python_requires='>=3.2',
    install_requires=[
        'biopython>=1.70',
        'colorama>=0.3.8',
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ),
)


