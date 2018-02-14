.. highlight:: shell

============
Installation
============


Stable release
--------------

To install Jade, run this command in your terminal:

.. code-block:: console

    $ pip install jade

This is the preferred method to install Jade, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for Jade can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/SchiefLab/jade

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/jadolfbr/jade/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ pip install -e . 

This will symlink your cloned code instead of copying it - very useful for development. 

Don't forget to source your bashrc/bash_profile or whatever other shell profile you have. 

Scripts
=======

All scripts and applications in jade/apps will be installed in your bin directory and available for use. 


.. _Github repo: https://github.com/SchiefLab/Jade
.. _tarball: https://github.com/SchiefLab/Jade/tarball/master
