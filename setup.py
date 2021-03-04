import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='stepRNA',
    version='1.0.2',
    author='Ben Murcott',
    author_email='bmm41@bath.ac.uk',
    description='Align short RNA seqeuncing reads to determine the length of of overhang.',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url='https://github.com/bmm514/stepRNA.git',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix"
    ],
    license='MIT',
    keywords=['Dicer', 'RNA', 'genomics', 'processing', 'sRNA', 'miRNA'],
    packages=["stepRNA"],
    python_requires=">=3.6",
    install_requires=[
        "pysam>=0.16.0.1",
        "biopython>=1.78",
        "numpy>=1.19.0"
        ],
    scripts=["bin/stepRNA",
    "stepRNA/stepRNA_run_bowtie.py",
    "stepRNA/stepRNA_cigar_process.py"
    ]
)

