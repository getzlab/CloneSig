from setuptools import setup, find_packages
import os

setup(
    name="CloneSig",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.0.0",
        "numpy>=1.19.0",
        "scipy>=1.5.0",
        "numpy-groupies>=0.9.0",
    ],
    python_requires=">=3.6",
    author="Liz Martin, Ignaty Leshchiner, Gad Getz",
    author_email="lmartin@broadinstitute.org",
    description="A tool for identifying convergent evolution across clones",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    url="https://github.com/getzlab/CloneSig",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    extras_require={
        'curveball': ['curveball @ git+https://github.com/getzlab/CurveBall.git'],
    },
) 