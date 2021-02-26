import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='stepRNA',
    version='0.1.11',
    author='Ben Murcott',
    author_email='bmm41@bath.ac.uk',
    description='Align short RNA seqeuncing reads to determine the length of of overhang.',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url='https://github.com/bmm514/stepRNA.git',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Licesnce :: OSI Approved :: MIT Licence",
        "Operating System :: Ubuntu",
    ],
    #packages=setuptools.find_packages(),
    license='MIT',
    packages=["stepRNA"],
    python_requires=">=3.6",
    install_requires=[
        "pysam>=0.16.0.1",
        "Bio>=0.3.0",
        "numpy>=1.20.1"
        ],
    scripts=["bin/stepRNA",
    "stepRNA/stepRNA-run_bowtie.py",
    "stepRNA/stepRNA-cigar_process.py"
    ]
)

