from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='jade',
      long_description=readme(),
      version='1.0',
      description='A repository for modules and applications to aid '
                  'in the design and analysis of Biological molecules, '
                  'especially when working with Rosetta or PyRosetta.',
      url='https://github.com/SchiefLab/Jade',
      author='Jared Adolf-Bryfogle',
      author_email='jadolfbr@gmail.com',
      license='BSD',
      packages=find_packages('jade'),
      install_requires=[
          'biopython',
          'numpy',
          'scipy',
          'scikit-learn',
          'pandas',
          'seaborn',
          'overrides'],
      classifiers=[
          'License :: OSI Approved :: BSD License',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Topic :: Utilities'],
      keywords='rosetta pyrosetta biology protein design bioinformatics carbohydrates pymol biochemistry modeling pdb',
      zip_safe=False,
      include_package_data=True,
      test_suite = 'nose.collector',
      tests_require = ['nose'],
)