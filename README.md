
Stochastic/deterministic invasion model 
===================================================

This is C++ code for spatially explicit simulation of an invasion process of a woody specie with a simple two stage structure juvenile/adults. 

The details of the models are in;

Saravia LA (2009) Algunas cuestiones sobre el espacio en ecología. Modelos, datos y aplicaciones Universidad Nacional de Buenos Aires (UBA).
<http://digital.bl.fcen.uba.ar/Download/Tesis/Tesis_4579_Saravia.pdf>


I am using GCC to compile it. For a crude and slow graphical output you will need the X11 SDL development libraries and the GRX graphics library http://grx.gnu.de/grx246um.htm

You also need the code from https://github.com/lsaravia/SpatialAnalysis.

### Files

SInvasions.cpp :    Deterministic and full stochastic model 
SInvasions01.cpp :    Deterministic + random error, linear juvenil mortality 
SInvasions02.cpp :    Deterministic + random error, Allee juvenil mortality 
SInvasions03.cpp :    Deterministic + random error, Allee juvenil mortality, inverse power dispersal
SInvasions04.cpp :    Deterministic + random error, Allee juvenil mortality, uniform long distance dispersal 

License
=======

	Copyright 2009 Leonardo A. Saravia
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
