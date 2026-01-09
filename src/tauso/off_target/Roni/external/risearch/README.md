# RIsearch
RIsearch: a tool for large-scale RNA–RNA, RNA-DNA, and DNA-DNA interaction prediction

## Installation

### RIsearch1
First, make sure that you have the following programs installed:

* gcc
* make

You can compile RIsearch1 by running the following from within the RIsearch1 folder.

	make RISEARCH

### RIsearch2

First, make sure that you have the following programs installed:

* gcc
* make
* CMake

Now, you should attempt to recompile the binary using the following command from
within the RIsearch2 folder:

	./rebuild.sh

If that does not work, and you have the gcc and CMake programs installed,
please report back to the author.

## Executing RIsearch
For more convenient use of RIsearch, add the installation folder or RIsearch1 and 
RIsearch2 to your PATH or copy the binary to a location that is in $PATH.

Here we provide examples on how to run RIsearch2 and RIsearch1 with default parameters.

See the manuals for RIsearch2 and RIsearch1 for the full list of options and examples.

### RIsearch2
First, generate the index structure for the target sequence(s) in the file target.fa
and store them in the file target.suf. 

	risearch2.x -c target.fa -o target.suf

To run RIsearch2 in default settings use the following command:

	risearch2.x -q query.fa -i target.suf

Note: the files query.fa and target.fa may contain several sequences, 
RIsearch will scan all vs. all

### RIsearch1
To run RIsearch1 in default settings use the following command:

	RIsearch -q query.fa -t target.fa

Alternatively, single sequences can be given directly on commandline with
-Q acgu -T acgu

Note: the files query.fa and target.fa may contain several sequences, 
RIsearch will scan all vs. all

## RIsearch manuals
Please read the full manuals for RIsearch1 and RIsearch2, included in the installation folders.

## Copyright

Copyright 2021 by the contributors (see RIsearch2/AUTHORS and RIsearch1/README files)

RIsearch1 and Risearch2 are released under the GNU General Public License
version 3. Note that libdivsufsort is packaged along with RIsearch.
Libdivsufsort comes with its own authors and copyright under
RIsearch2/libdivsufsort-2.0.1/{AUTHORS,COPYING}

GNU GENERAL PUBLIC LICENSE


This is a free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License, either version 3 of the License, or
(at your option) any later version. You should have received a copy of the GNU General Public License
along with RIsearch, see file COPYING. If not, see <http://www.gnu.org/licenses/>.

This software is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citations

If you use RIsearch in your publication please cite:

**RIsearch2: suffix array-based large-scale prediction of RNA-RNA interactions and siRNA off-targets**
Alkan F, Wenzel A, Palasca O, Kerpedjiev P, Rudebeck AF, Stadler PF, Hofacker IL, Gorodkin J *Nucleic Acids Res*. 2017 May 5;45(8):e60

**RIsearch: fast RNA-RNA interaction search using a simplified nearest-neighbor energy model**
Wenzel A, Akbasli E, Gorodkin J. *Bioinformatics*. 2012 Nov 1;28(21):2738-46. Epub 2012 Aug 24

## Contact

In case of problems or bug reports, please contact: <software+crispron@rth.dk>

