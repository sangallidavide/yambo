## The Ydriver Yambo Library

As part of the modularization of the Yambo code, one of the goals is to split the source into different blocks (or modules), and eventually promote them into libraries. So far we created, as an experimental library, the driver library.
The source code of the library is not contained anymore in the main repo of the yambo code, but it has been moved in this dedicated repository.

### Yambo
Yambo doesn't stand for anything like "Yet Another Many-Body cOde", for instance.  Unless you really want it to. 

For more information, see the general documentation on the web site http://www.yambo-code.eu/

Compiling:  http://www.yambo-code.eu/download/compiling.php
Installing: http://www.yambo-code.eu/download/installing.php
Getting started: http://www.yambo-code.eu/doc/getting_started.php

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
for more details (see COPYING file).

### Specific compilation flags

Nvidia 
./configure FC=nvfortran FPP="cpp -traditional" FCMFLAG="-Mnomain"
