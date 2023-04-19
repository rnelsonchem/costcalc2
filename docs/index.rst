*costcalc2* 
===========
A Chemical Synthetic Route Cost Calculation Package
---------------------------------------------------

Welcome to the *costcalc2* documentation. This Python3 package provides tools
for calculating the overall raw material costs (RMCs) and PMI for chemical
synthetic routes of arbitrary complexity. 

There is a `companion web application <https://costcalc.rnelsonchem.com/>`_
that can be used to run these cost calculations without writing any Python
code. (Note: because of how the app is currently hosted, the web app may be
slow to load.) Because the web application is the likely first entry point for
new users looking to use this tool, the first section below will provide
details on creating properly formatted materials and reactions tables in
Excel, which can then be uploaded and processed in the web app. The second
section will describe how to interpret the output table from these
calculations.

Subsequent sections will be geared towards users that are interested in coding
the cost calculations directly in Python. This provides more flexibility
in the specifics of the user input/output and access to additional features
that are not implemented in the web application.

Installation
------------

The *costcalc2* package can be installed using either ``pip``

.. code-block:: console 

    $ pip install costcalc2
    
or ``conda``

.. code-block:: console 

    $ conda install -c rnelsonchem costcalc2


Contents
--------

.. toctree::
   :maxdepth: 2

   tables_basics
   interpret_results

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
