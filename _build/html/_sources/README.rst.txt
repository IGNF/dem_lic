===========================
DEM Generalization with LIC
===========================

A Python package for Digital Elevation Model (DEM) generalization, leveraging Line Integral Convolution (LIC) and adaptive Gaussian blur techniques to simplify terrain representations while preserving key features. The method processes large DEM files efficiently by dividing them into overlapping blocks.

Documentation
-------------

Full documentation is available at:

.. code-block:: text

    [Link to Sphinx Documentation - Coming Soon]

What does it do?
----------------
This algorithm generalizes Digital Elevation Models (DEMs) by combining adaptive Gaussian blur based on terrain curvature, flat and steep area detection for feature enhancement, and Line Integral Convolution (LIC) to emphasize ridges and key terrain features. Iterative processing ensures improved visual coherence and clarity.

.. image:: docs/images/dem_to_generalization_2m.png
   :alt: Exemple of generalization


Installation
------------

To install the package, clone the repository and install it in editable mode.

HTTP method ::

    git clone https://github.com/ESaint-Denis/dem_lic.git
    cd dem_lic
    pip install -e .

SSH method ::

    git clone git@github.com:ESaint-Denis/dem_lic.git
    cd dem_lic
    pip install -e .

Command Line Quick Start
------------------------

The command-line interface (CLI) for **dem_lic** provides an easy way to process Digital Elevation Models (DEMs) using the generalization algorithm. The CLI supports input and output file specifications, as well as customization of processing parameters.

.. code-block:: bash

    $ dem_lic input_dem.tif output_generalized_dem.tif

Optional arguments can be used to adjust the process:

.. code-block:: bash

    $ dem_lic input_dem.tif output_generalized_dem.tif --n_iterations 5 --overlap 20

For detailed usage instructions and a complete list of options, refer to the `CLI Documentation <http://example.com/cli-docs>`_.

Python Quick Start
------------------

You can also use **dem_lic** directly within your Python scripts. Here is a quick example:

.. code-block:: python

    from dem_lic.generalization import generalization

    generalization(
        MNT_input_path="path_to_input_dem.tif",
        output_path="path_to_output_dem.tif",
        block_size=2000,
        overlap=20,
    )

License
-------

This project is licensed under the MIT License. See the LICENSE.TXT file for details.


