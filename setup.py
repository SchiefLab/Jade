from setuptools import setup, find_packages
import os,glob

if hasattr(__builtins__, 'raw_input'): input = raw_input

def readme():
    with open('README.md') as f:
        return f.read()

def get_all_scripts_to_install(public_dir='apps/public', pilot_dir='apps/pilot'):

    #Public Apps

    possible_apps = ['\tall', 'none']
    apps_to_install = []
    all_apps =(sorted([os.path.basename(d) for d in glob.glob(os.path.join(public_dir, "*")) if os.path.isdir(d)]))

    possible_apps.extend(all_apps)

    app_modules = input("\nPlease enter any apps to install (space separated): \n\n"+"\n\t".join(possible_apps) +"\n\n")
    app_modules = str(app_modules).strip().split()

    if app_modules.count('all'):
        app_modules = all_apps
        app_modules = [public_dir + "/" + d for d in app_modules]

    elif app_modules.count('none'):
        app_modules = []
    else:
        app_modules = [public_dir + "/" +d for d in app_modules]


    #Pilot Apps
    possible_pilot_apps = ['\tall', 'none']
    pilot_apps = (sorted([ os.path.basename(d) for d in glob.glob(os.path.join(pilot_dir, "*")) if os.path.isdir(d)]))
    possible_pilot_apps.extend(pilot_apps)

    print(pilot_apps)

    pilot_app_modules = input("\nPlease enter any pilot apps to install (space separated):\n\n" + "\n\t".join(possible_pilot_apps)+"\n\n")
    pilot_app_modules = str(pilot_app_modules).strip().split()

    if pilot_app_modules.count('all'):
        pilot_app_modules = pilot_apps
        pilot_app_modules = [pilot_dir + "/"+d for d in pilot_app_modules]
    elif pilot_app_modules.count('none'):
        pilot_app_modules = []
    else:
        pilot_app_modules = [pilot_dir + "/"+d for d in pilot_app_modules]

    #Get all app code paths.
    app_modules.extend(pilot_app_modules)

    all_app_paths = []
    for app_module in app_modules:
        apps = glob.glob(os.path.join(app_module, "*.py"))
        all_app_paths.extend(apps)
    return all_app_paths

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
      scripts = get_all_scripts_to_install(),
)