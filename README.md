# QTM Generator

This script makes a Quaternary Triangular Mesh (QTM) to tessellate the planet geodetically into a discrete global grid system (DGGS) based on an octohedron. Each triangular facet resursively subdivides into four smaller ones to build a hierarchical mesh. We directly implement the geometry and tile indexing scheme developed by Geoffrey Dutton, and described in the following publication:

Dutton, Geoffrey H. "Planetary Modelling via Hierarchical Tessellation." In Procedings of the AutoCarto 9 Conference, 462–71. Baltimore, MD, 1989. https://pdfs.semanticscholar.org/875e/12ce948012b3eced58c0f1470dba51ef87ef.pdf



Written by Paulo Raposo (pauloj.raposo [at] outlook.com) and Randall Brown (ranbrown8448 [at] gmail.com) at the [Department of Geography, University of Tennessee, Knoxville](http://geography.utk.edu/).

## Usage:

All geometry is defined in [WGS 84](http://spatialreference.org/ref/epsg/4326/). The script outputs to a GeoJSON file for each hierarchical level in the QTM, making as many levels as the user requests.

The output folder specified must pre-exist. Output files are named `qtmlvlX.geojson`, where `X` is the integer number of the level, 0 or higher.

### Examples on the command line:
Unix-like: `python qtmgenerator.py /home/username/qtmfolder 12`

## Dependencies: 

* nvector (see https://pypi.python.org/pypi/nvector and http://www.navlab.net/nvector), 
* OGR Python bindings (packaged with GDAL, see https://pypi.python.org/pypi/GDAL).
