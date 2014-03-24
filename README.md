# phyloconverge

## DESCRIPTION

`phyloconverge` has not been tested extensively and is being actively developed. Please use with caution. I appreciate hearing about your experience with the program.

`phyloconverge` is an algorithm to look for convergence based on a particular trait in a phylogenetic tree. The traits can be quantitative (e.g., size in centimeters or depth of observation) or can be categorical (e.g., associated with hydrothermal vents or not associated with hydrothermal vents).

## AVAILABILITY

https://github.com/josephryan/phyloconverge (click the "Download ZIP" button at the bottom of the right column).

### DEPENDENCIES

NOTE: python is called within perl to process the treefile

General system tools:
- [Perl] (http://www.cpan.org/), which comes with most operating systems
- [Python] (http://www.python.org/)

Additional libraries
- [DendroPy Phylogenetic Computing Library] (http://pythonhosted.org/DendroPy)

## INSTALLATION

To install `phyloconverge` and documentation, type the following:

    perl Makefile.PL
    make
    make test
    sudo make install

## RUN

    phyloconverge --treefile=NEWICKTREEFILE --tabfile=TABFILE_W_ATTRIBUTES

## EXAMPLE ANALYSES
    
    phyloconverge --treefile=examples/converged.tre \
                  --tabfile=examples/depths.txt

    phyloconverge --treefile=examples/notconverged.tre \
                  --tabfile=examples/depths.txt

## DOCUMENTATION

Extensive documentation is embedded inside of `phyloconverge` in POD format and
can be viewed by running any of the following:

        phyloconverge --help
        perldoc phyloconverge
        man phyloconverge  # available after installation

## COPYRIGHT AND LICENSE

Copyright (C) 2013 Joseph F. Ryan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program in the file LICENSE.  If not, see
http://www.gnu.org/licenses/.
