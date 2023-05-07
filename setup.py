# coding utf8
import setuptools
from snpfilter.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="SNPfilter",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="A handy little tool for filtering SNPs.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/SNPfilter",

    entry_points={
        "console_scripts": ["SNPfilter = snpfilter.cli:main"]
    },  

    packages=setuptools.find_packages(),

    install_requires=[
        "pysam>=0.21.0",
        "toolbiox>=0.0.39",
    ],

    include_package_data=True,

    python_requires='>=3.5',

)