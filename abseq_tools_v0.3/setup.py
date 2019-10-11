import setuptools


if __name__ == '__main__':
    setuptools.setup(
        name='bdgutils',
        version='0.3',
        packages=setuptools.find_packages("src"),
        package_dir={'': 'src'},
        package_data={'':['AbSeqV2Final.fasta', 'AbSeqV2Final.gtf','BDSampleTags.fasta','BDSampleTags.gtf',
                          'BDSampleTagsMM.fasta','BDSampleTagsMM.gtf']},
        entry_points={
            'console_scripts': [
                'abseq_counts = bdgutils.abseq_counts:main',
                'demux_gene_counts = bdgutils.demux_gene_counts:main',
                'create_ref = bdgutils.create_ref:main'
            ],
        },
        install_requires=['pandas', 'numpy', 'scipy', 'begins', 'tables', 'pysam', 'matplotlib'],
        setup_requires=[],
        tests_require=['pytest'],
        author='Janice Lai',
        author_email='janice.h.lai@bd.com',
        description='Utilities for processing data from BD AbSeq',
        long_description="This package contains tools for processing data generated using the BD AbSeq",
        license = 'CC BY-NC-SA'
)
