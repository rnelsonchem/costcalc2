from setuptools import setup, find_packages

with open('README.md') as file:
    long_description = file.read()

setup(
    name = "costcalc2",
    version = "0.3",

    description = "Cost calculations for synthetic chemical preparation routes",
    url = "https://github.com/rnelsonchem/costcalc2",
    long_description = long_description,

    author = "Ryan Nelson",
    author_email = "rnelsonchem@gmail.com",

    license = "GNUv3",
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU License',
        'Programming Language :: Python :: 3',
    ],

    keywords = "Chemical Synthesis Cost Pricing",

    packages = find_packages(),
    install_requires = [
        'numpy>=1.17',
        'pandas>=0.25',
        'matplotlib>=2.0',
        'openpyxl>=3.0',
        ]

)


