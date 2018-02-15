
from setuptools import setup, find_packages
import os,glob,re



def readme():
    try:
        with open('README.rst') as f:
            return f.read()
    except IOError:
        return ''

def get_all_scripts_to_install_user(public_dir='apps/public', pilot_dir='apps/pilot'):
    """
    An attempt at user-level control of the install.  Does not work.  Pip hangs.
    The function itself has been tested outside of this.
    :param public_dir: 
    :param pilot_dir: 
    :return: 
    """
    if hasattr(__builtins__, 'raw_input'): input = raw_input

    print "Gathering scripts to install...\n"

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

def get_all_scripts_to_install(public_dir='apps/public', pilot_dir='apps/pilot'):
    all_scripts = []
    for outer in [public_dir, pilot_dir]:
        for app_dir in (sorted([ d for d in glob.glob(os.path.join(outer, "*")) if os.path.isdir(d)])):
            #print "reading " + app_dir
            f = glob.glob(app_dir+"/"+"*.py")
            #print(f)
            for script in f:
                if not re.search("__init__", script):
                    all_scripts.append(script)

    #print all_scripts
    print "Found scripts:"
    print "\n\t".join(all_scripts)
    print "\n\n"
    return all_scripts

def find_all_packages():
    #print "Finding Packages."
    p = ['jade/'+ sub for sub in find_packages('jade')]
    #print repr(p)
    return p

setup(name='bio-jade',
      long_description=readme(),
      version='1.0',
      description='A repository for modules and applications to aid '
                  'in the design and analysis of Biological molecules, '
                  'especially when working with Rosetta or PyRosetta.',
      url='https://github.com/SchiefLab/Jade',
      author='Jared Adolf-Bryfogle',
      author_email='jadolfbr@gmail.com',
      license='BSD',
      packages= find_all_packages(),
      install_requires=[
          'biopython',
          "sphinx-argparse",
          'numpy',
          'scikit-learn',
          'pandas',
          'seaborn',
          'overrides',
          'weblogo',
          'prettytable'],
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
      test_suite='tests',
      scripts = get_all_scripts_to_install(),
)