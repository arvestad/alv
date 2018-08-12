import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('alv/version.py') as fh:
    exec(fh.read())
    
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


