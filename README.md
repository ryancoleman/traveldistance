![Travel Depth](http://crystal.med.upenn.edu/tdimages/traveldepth.png)

[![doi](https://zenodo.org/badge/3853/ryancoleman/traveldistance.png)](http://dx.doi.org/10.5281/zenodo.10213)


traveldistance
==============

Github version of Travel Distance code by me (Ryan Coleman) and my advisor (Kim Sharp). 
See http://crystal.med.upenn.edu/travel_distance.html for official download, academic paper references, license, etc.

INSTALLING
==========

This is a guide to installing Travel Depth, CHUNNEL, Burial Depth, and CLIPPERS
 and the various dependencies it has. Please remember to cite the 
following if you use this work in your research:

Travel Depth:

Coleman, RG, Sharp, KA. Travel Depth, A New Shape Descriptior for 
Macromolecules: Application to Ligand Binding. Journal of Molecular 
Biology 362(3), pp. 441-458, 22 September 2006. 
http://dx.doi.org/10.1016/j.jmb.2006.07.022

CHUNNEL:

Coleman, RG, Sharp, KA. Finding and characterizing tunnels in 
macromolecules with application to ion channels and pores. Biophysical 
Journal. 96(2). pp 632-645. January 2009. 
http://dx.doi.org/10.1529/biophysj.108.135970

Burial Depth:

Coleman, RG, Sharp, KA. Shape and evolution of thermostable protein 
structure. Proteins: Structure, Function, and Bioinformatics. 78(2). pp 
420-433. 1 February 2010. 
http://dx.doi.org/10.1002/prot.22558

CLIPPERS:

Coleman, RG, Sharp, KA. Protein Pockets: Inventory, Shape, and Comparison. 
J Chem Inf Model. 50(4). pp 589-603. 26 April 2010.
http://dx.doi.org/10.1021/ci900397t

For additional information, or the original source of these files:

http://crystal.med.upenn.edu
http://ryancoleman.name

0. Systems

All code is developed and tested on GNU/Linux systems, and has been tested 
on recent Mac OS X systems. It may work under cygwin on Windows systems.

1. External Downloads/Installs

The Travel Depth Code relies on Python 2.5 (it may work with 2.6 but is 
untested, a port to the 3 branch is planned). Also, visualization support is 
provided by a python script for PyMOL. CLIPPERS requires numpy or Numeric, 
it is optional for Travel Depth, Burial Depth and CHUNNEL. 

These 3 software packages may be installed on your system already, or if
necessary, download and install them here:

http://www.python.org/download/releases/2.5.2/

http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=1351

http://pymol.sourceforge.net/

You must set the $TDHOME environmental variable to this directory (the 
directory with bin/ and src/ among others.

The included python scripts run using the python2.5 command. You must have 
this set and in your path unless you want to call each script by saying 

`/path/to/python2.5 $TDHOME/src/tstTravelDepth.py file.pdb`

for instance. Again, if you've set python2.5 to call the correct 
interpreter, you can just run the scripts like this:

`$TDHOME/src/tstTravelDepth.py file.pdb`

Before running these scripts, you must also install the required software.

2. Quickhull

Quickhull or qhull is redistributed here, available at
http://www.qhull.org as well, and it is required for the Travel Depth
code to run. If using the supplied version, install as follows:

`cd $TDHOME/qhull/src`
`make install`

$TDHOME/bin/qhull is the qhull executable, other qhull executables are also
placed there. 

An alternative code khull is under development and included since qhull 
does not work on some 64 bit linux installations. khull is compiled in the 
next step, after completing the compilation you can rename khull to qhull 
in the $TDHOME/bin directory. Alternatively, you can change the flag in 
tstConvexHull.py on line 12 to runKHull=True.

3. Surface Generation code (TST)

Travel Depth can work with any triangulated surface, but currently
support only for surfaces generated with the included code and the 'tst'
file format. Install as follows:

`cd $TDHOME/tst`
`make all`

Which compiles and puts the executables into ./bin 

4. If you want to analyze CLIPPERS data, you will most likely need to 
install aiSee. aiSee is free for academic use and can be downloaded from:

http://www.aisee.com/download/

For almost any operating system. aiSee imports .gdl files, which are 
created by CLIPPERS as described in the RUNNING_TRAVEL_DISTANCE.txt file.

5. All python scripts are are in pymol/ or src/, and need no
installation, unless you want to add the directories to your path or run
the tstdisp.py script everytime PyMOL starts.

There are several scripts to run, described in the 
RUNNING_TRAVEL_DISTANCE.txt file

RUNNING
=======

This walks through running an example of the scripts you can run.
Ryan G. Coleman, Kim Sharp Lab, http://crystal.med.upenn.edu

Please remember to cite the following if you use this work in your research:

Travel Depth:

Coleman, RG, Sharp, KA. Travel Depth, A New Shape Descriptior for
Macromolecules: Application to Ligand Binding. Journal of Molecular
Biology 362(3), pp. 441-458, 22 September 2006.
http://dx.doi.org/10.1016/j.jmb.2006.07.022

CHUNNEL:

Coleman, RG, Sharp, KA. Finding and characterizing tunnels in
macromolecules with application to ion channels and pores. Biophysical
Journal. 96(2). pp 632-645. January 2009.
http://dx.doi.org/10.1529/biophysj.108.135970

Burial Depth:

Coleman, RG, Sharp, KA. Shape and evolution of thermostable protein
structure. Proteins: Structure, Function, and Bioinformatics. 78(2). pp
420-433. 1 February 2010.
http://dx.doi.org/10.1002/prot.22558

CLIPPERS:

Coleman, RG, Sharp, KA. Protein Pockets: Inventory, Shape, and Comparison.
J Chem Inf Model. 50(4). pp 589-603. 26 April 2010.
http://dx.doi.org/10.1021/ci900397t

For additional information, or the original source of these files:

http://crystal.med.upenn.edu
http://ryancoleman.name

RUNNING TRAVEL DISTANCE:

1. As described in installing, python2.5 must be in your path and must be 
the python2.5 executable with numpy or Numeric installed. $TDHOME must be 
set to the directory this file is in, since $TDHOME/bin and $TDHOME/src 
are referenced by several files.

2. The input is a pdb file, you should remove heteroatoms (depending on
what your desired analysis is like). Waters should be ignored but you may as
well remove them too. From here on out, I'll refer to pdb file as file.pdb

3. It is best to work in a new directory with your pdb file in it, there 
are a lot of temporary files generated at various points. These should be 
properly cleaned up, but it is best to not clutter things. Note that 
running many computations at once must be done in separate directories (on 
multi-core machines for instance) due to temporary file writing by the tst 
suite of programs.

4. tstTravelDepth.py is the script that creates tst and phi files and 
actually calculates travel depth. Run it as 
`$TDHOME/src/tstTravelDepth.py file.pdb`. The additional commandline 
options are [gridsize] [tri|mesh] [probe] [radius scale] [pathToExecs], 
which must be in that order. [gridsize] is the size in angstroms per grid 
side, tri or mesh denote 2 different surface generation techniques (tri is 
a gaussian emulation of probe sphere rolling, mesh is actually rolling 
probe spheres). [probe] is the probe size, and should not be changed 
outside of reasonable ranges of water probes, and is not completely 
accurate when using the tri meshing scheme. [radius scale] can change the 
van der Waals radii but a constant multiplier, for instance 0.9 would 
cause all the atoms to be 10% smaller. [pathToExecs] should not need to be 
used but if the convex hull and fortran code are not in your path you will 
need to supply the $TDHOME/bin directory here.

If everything works, you should get a file.nocav.tst and file.nocav.phi 
file, and the travel depth results are in the .tst file.

5. This step is to visualize the travel depth results. In PyMOL, at the 
prompt, `run $TDHOME/pymol/tstdisp.py` then the functions tstOpen and tstDraw 
are available. `tstOpen file.tst` will load a file then `tstDraw` 
will show the standard depth coloring on the surface. The next most useful 
option is `tstDraw depth,maxval=XX` which changes what blue corresponds to 
in travel depth, most useful for coloring several figures according to the 
same color scale. Also try `tstDraw depthmono` or `tstDraw depthmonorev` 
if you'd like black and white figures.

Additionally there is a new function tstOpenDraw that opens any number of 
files by calling it as `tstOpenDraw filemask, depth` or other options 
after depth. calling `tstOpenDraw *.tst, depth` will produce surfaces for 
all tst files in the current directory.

Other options are documented in the tstdisp.py file itself.

6. This step walks you through a CHUNNEL run, assuming you already have 
file.nocav.tst and file.nocav.phi files produced by running travel depth 
on your file.pdb of interest. 

An optional step is to run the following `$TDHOME/src/tstTravelDist.py 
countholes file.nocav.tst file.nocav.phi` which simply counts the number 
of holes in the surface according to the Euler characteristic. If there 
are no holes, running CHUNNEL will produce no results, if there are lots 
and lots of holes, running CHUNNEL will be very slow. In either case, try 
adjusting the probe radius when constructing the surface.

Run `$TDHOME/src/tstChunnel.py file.nocav.tst file.nocav.phi` to 
start CHUNNEL. For extremely large or hole-y proteins, this may take quite 
awhile, it is recommended to run in the background and log the results to 
a file like this: `nohup $TDHOME/src/tstChunnel.py file.nocav.tst 
file.nocav.phi >& chunnel.log.txt &`

The output log file file.nocav.tst.findholes.log is a space-delimited file 
containing data on the paths found. It can be easily imported into Excel 
or your favorite spreadsheet software. The following is a list of the 
columns "number endNumOne endNumTwo plugs stepLength pathLength 
pathMinRadius pathMaxInsi deRadius endMinOne endMinTwo minimaCount 
travelDepthMax windingMetric avgTheta"

Obviously for a more complete understanding of the following brief 
descriptions, you'll need to read the CHUNNEL paper.

number is the reference number of each hole, starting at 1. This is used 
to write the additional files (described later.)

endNumOne and endNumTwo are the numbers of each end of the path. These are 
arbitrary but consistent. If two paths share an end number they share 
an endpoint.

plugs is a list of the 'plugs' the path encounters.

stepLength is the number of grid points along the path

pathLength is the distance the path travels

pathMinRadius is the minimum radius along the path

pathMaxInside is the maximum between the first minima from either end.

endMinOne and Two are the lowest and second lowest minima

minimaCount is the number of unique minima in terms of radius along the path

travelDepthMax is the maximum travel depth along the path

windingMetric is the pathLength/distance(start, end) of the path. higher 
is windier

avgTheta is the average theta between adjacent path points. Not very 
useful (not used in paper).

Output Files:

file.nocav.tst.outside.py is a pymol file that contains spheres at each 
entry/exit point

file.nocav.tst.tree.py is a pymol file that contains spheres at each grid 
point in the entire tree of paths

And for each hole number, like hole 1:

file.nocav.tst.1.pore.py is a pymol file with the radius shown

file.nocav.tst.1.path.py is a pymol file where each grid point on the path 
is just a small dot

file.nocav.tst.1.radii.txt is just a list of the radii

file.nocav.tst.1.residues.pathmin.pdb is a pdb file of the residues near 
the choke point

file.nocav.tst.1.residues.path.pdb is a pdb file of the residues lining 
the entire path

7. This step walks you through a Burial Depth run. This requires running 
Travel Depth, but not CHUNNEL. So you already have a file.nocav.tst and 
file.nocav.phi set of files.

Run burial depth as `$TDHOME/src/tstBurialDepth.py file.nocav.tst 
file.nocav.phi`. It will append data to the tst file but most importantly 
it will create a file.nocav.tst.mesh.travelin.pdb file that has the burial 
depth in Angstroms replaced in the b-factor column. 

You can analyze this data in some way (importing the pdb file into excel, 
averaging by residue, etc) or you can visualize it by importing it into 
pymol and choosing color -> spectrum -> b-factors. Other visualization 
programs may also have ways to visualize b-factor data on a color scale.

8. This sectino describes running CLIPPERS. Again, this only depends on 
running Travel Depth, so you have a file.nocav.tst and file.nocav.phi file 
already created.

For every file.nocav.tst and file.nocav.phi file you have, you want to map 
the pockets of the surface. This is done by running 
`$TDHOME/src/tstPocketMap.py file.nocav.tst file.nocav.phi` Again, like 
CHUNNEL, this can take awhile for large proteins, so you may want to run 
it like this: `nohup $TDHOME/src/tstPocketMap.py file.nocav.tst 
file.nocav.phi >& file.nocav.pocketmap.log.txt &` 

After this has completed, the file.nocav.tst has been modified and a new 
file called file.nocav.tst.tree.tm3 has been created. First, to visualize 
the pockets in PyMOL.

In PyMOL run $TDHOME/pymol/tstdisp.py which makes new commands available.

Next, in PyMOL, run `tstOpen file.nocav.tst` 

Now, individual pockets can be displayed by running `tstDraw group=XXX` 
where XXX is a specific group number. The file.nocav.tst.tree.tm3 has a 
list of pockets, the first column is the pocket or group ID number. 

Alternatively, you can run the command `tstDrawDeepest` which draws the 
deepest pocket found. If you run `tstDrawParent` after this command or 
after drawing any given pocket, the parent of this pocket (up the tree, 
bigger, less deep) will be displayed. If you run `tstDrawParent 3`, the 
great-grandparent pocket will be displayed, replace the 3 to go even 
further up the tree at once. You cannot go down the tree easily since it 
is branched but these commands allow you some control over navigating the 
tree structure of the pockets. 

Alternatively, explore the file.nocav.tst.tree.tm3 file to find a pocket 
of interest and display it using the `tstDraw group=XXX` command. Future 
work will include utilities to pick out interesting pockets for you or 
label pockets as bottlenecks or clefts.

As described in the CLIPPERS paper, pockets can be compared or clustered 
in many ways. One way to do this is as follows, assuming you have 2 tm3 
files file1.tm3 and file2.tm3:

`$TDHOME/src/tm3search.py file1.tm3 file2.tm3`

The output from this will result in a series of .gdl files at various 
numbers of connections added to the clustering. These are viewable using 
aiSee. Simply go to File->Open within aiSee and select a .gdl file of 
interest. Note the black bars (representing nodes that have similar 
shapes) and grey arrows (parent->child links). Nodes are colored according 
to the adjacent nodes' conservation, with red being conserved and blue 
being unconserved. If you use the 'information' tool in aiSee and select a 
node you will get a dialog box containing a string that looks like this

`tstOpenDraw file.nocav.tst, group=45` which you can copy and paste into 
PyMOL to display that particular node.

If you'd like more control over the parameters used when comparing whole 
trees of pockets, you can use the following script 

`$TDHOME/src/tm3searchParameters.py alpha minVolume maxVolume 
startConnections endConnections stepConnections 
useNormalMeansStds(boolean) [tm3 file list]`

The parameters are described fully in the CLIPPERS paper, but briefly:

alpha is the parameter used to control the redundancy of pockets in the output

minVolume is the parameter of the smallest volume in cubic angstroms that 
a pocket must be to be considered

maxVolume is the maximum volume that will be considered

The connections parameters define how many connections will be added to 
the cluster before the series of .gdl files for aiSee are created.

The last parameter is whether to use the normal means and standard 
deviaitions from a 'standard' set of proteins when computing z-scores or 
whether to calculate the means and standard deviations from the actual 
pockets in the tm3 files used as input.

Finally is the list of tm3 files. Any number of tm3 files can be used, but 
this is very ram-limited as all files and connections between pockets are 
stored in memory in this implementation. 20 tm3 files with about a 
thousand pockets each uses about 2 gigabytes of RAM and is the limiting 
factor, such a large run can take about a day. Comparing 2 such files 
takes a few minutes.

Another option (fully documented in the CLIPPERS paper) is to find the 
pockets that best overlap a set of residues, then compare their geometric 
shapes across a series of structures. This is done with the following command:

`$TDHOME/src/tm3rescomp.py [name|number] residueList tm3List`

Either "name" or "number" should be specified depending on whether the 
list of residue ids is a named and numbered list "ALA17+ARG18+CYS20+LEU30" 
or just a numbered list "17+18+20+30". The residue list should be 
specified next in + delimited form and finally any number of tm3 files can 
follow. The output is the pocket numbers and snippets of code to display 
them, then a matrix of the distances between each pocket found.

Running this same command as the following:

`$TDHOME/src/tm3resRefineComp.py [name|number] residueList tm3List`

Runs the refinement algorithm, where pockets are found that are the 
minimum shape distance from each other as well as having at least 50% 
overlap with the set of residues. Again, this is fully described in the 
paper and the output is the same.

The command to cluster and examine the clusters (as is done to identify 
enzyme shapes) is the following:

`$TDHOME/src/tm3searchAlphaExamineClusters.py alpha tm3fileList` 

which is similar to the tm3search.py except it outputs data about the 
clusters.

The final included command is tm3summarize.py:

`$TDHOME/src/tm3summarize.py tm3fileList`

Which summarizes the volume, surface area and number of mouths across the 
tm3 files given. This file can be easily modified to include data on any 
geometric property of interest.
