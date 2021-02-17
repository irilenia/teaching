
Author: Irilenia Nobeli

Date: 17/02/2021

# Chemoinformatics with OpenBabel
## A short introduction to OpenBabel
OpenBabel is an open source chemical toolbox written in C++ designed for searching, converting, analysing and storing chemical data. Before you do anything else, I suggest you find out a bit more about OpenBabel using the [wiki](http://openbabel.org/wiki/Main_Page) pages dedicated to the project. In particular, take a look at the [About OpenBabel](http://openbabel.org/wiki/Open_Babel:About) pages.

OpenBabel can be used in different ways:
* it comprises [ready-to-use programs](http://openbabel.org/wiki/Guides) which can be run on their own without you needing to write extra code
* you can use the OpenBabel library of classes in your own scripts, e.g. you can write Perl (or Python) scripts       that will use OpenBabel's programmatic interface via a Perl wrapper.
* you can write your own C++ programs and use the OpenBabel library of classes directly

Before you proceed further, it's a good idea to familiarise yourselves a bit with the types of classes available and what OpenBabel is capable of. Take a quick look at the relevant pages on the wiki ([capabilities](http://openbabel.org/wiki/Capabilities) and [classes/API](http://openbabel.org/api/)). Using the OpenBabel classes, however, is outside the scope of this tutorial.

### How to access OpenBabel
**Ignore the paragraph in italics if you are connecting from home.** Instead please follow the guidelines I posted on Moodle to connect to thoth.

*Please read the following carefully. This part of the practical DOES NOT NEED YOU to install anything yourselves. However, if you want to know how to install the OpenBabel library, you can consult the installation guides online. OpenBabel is not installed on the College machines. Instead you need to be able to ssh via a terminal to a server in the department. For this practical you should connect to the server ssh-ext.cryst.bbk.ac.uk . To login to a departmental server from the College computers, you need to use the same procedure you used when you did your Python tutorials. If this process doesn't work, you can use __PuTTY__ but this is occasionally problematic when the commands are really long and go beyond a single line, so it's best avoided for this tutorial. Remember that to login on "pandora" you will need your departmental login and password, not the College (ITS) ones.*

Once you are logged on to the correct departmental server and you are at the directory you want to work in (this could be your home directory or a directory provided to you by the department on a separate network drive (see my Moodle guidelines), create a directory where you can put your files relating to this practical, e.g.
```
mkdir chemoinfo_practical
```
The ready-to-use programs of OpenBabel are installed under /usr/bin and so running them should be straightforward. You should not need to include the whole path (/usr/bin/babel) as the path to /usr/bin must already be included in your PATH variable. To check that you can indeed call these binaries try:
```
babel -H
```
If this does not produce a help screen for babel then please let us know and we'll look into it. Now you are ready to do the practical.

## Ready-to-use openbabel programs - babel and obgen
Before using any program, try to find out as much as possible about it from the OpenBabel wiki pages. Then run each program from the command line (using the executables installed in the directory I told you about in the Introduction page).

We'll start by using the program babel to perform  __chemical file translations__. If you want to visualise the results you have various options. Here are a few:
1. If you have a SMILES string you can use [Molview](http://molview.org) or [Molinspiration's Galaxy](http://molinspiration.com/cgi-bin/galaxy) structure generator facility.
2. You can download (at home) the [CACTVS](http://www.xemistry.com/) toolkit and use the editor csed to view most of the available chemical data files (again no need for coordinates to be provided). This is possibly the most reliable editor that is available free to academics.
3. Alternatively you can use the web sketcher available at: [Xemistry](http://85.214.192.197/edit/frame.html).
4. If you have 2D or 3D coordinates for your molecule, you can use __Chimera__ to view some of the more common chemical structure files.

In this practical we will work with the phosphodiesterase 5 inhibitor __Viagra__. Create a file called viagra.smi (eg using an editor like nano, nedit or vi) and put in it the following SMILES string:
```
CCCc1nn(C)c2C(=O)NC(=Nc12)c3cc(ccc3OCC)[S](=O)(=O)N4CCN(C)CC4
```
Then use the program __babel__ to translate this into an SDF file called __viagra.sdf__. Open the resulting file with a text editor and check the coordinates section of the atoms. _What do you notice_?

Try the same translation from SMILES to SDF but this time add :
1. a title (use for example "Viagra") and
2. hydrogens at physiological pH (try 7.4).

_What is the charge of the molecule produced?_ (hint: there should be a CHG line in the SDF file) Think about how this charge has arisen and whether you can trust the result.
Visualise all the above results using one of the options I listed above or any other facility you wish to use on the web or locally.

As you started with a SMILES string you will not have any coordinates for your molecule. The program __obgen__ allows you to produce optimised 3D coordinates (note that the obgen program produces the final file in the standard output and so you need to redirect the output to a file using __>__). Do this for the SMILES string above (this might take a few minutes to execute as multiple conformations need to be produced and scored). Visualise the result in __Chimera__ so you can check that a reasonable 3D geometry has been obtained (if you are working from home, you will need to transfer the file to your home computer using scp or an application such as __FileZilla__; if you are working on a College machine, you will need to copy the file from the departmental server to the College machines using any file transfer/ftp facility provided by the College e.g. __FileZilla__).

## Ready-to-use openbabel programs - Searching for substructures
For this part of the tutorial, we will use as input an MDL mol/sdf formatted file (one of the most commonly used formats) which can be downloaded from here: [1equ_ligonly_edited.sdf](http://people.cryst.bbk.ac.uk/~ubcg71a/OLD/teaching/chemoinfo_practical2b/1equ_ligonly_edited.sdf). Alternatively on _thoth_ or other departmental machine you can copy this file from:
```
/d/in4/u/ubcg71a/teaching/chemoinfo/1equ_ligonly_edited.sdf
```

This file contains the results (coordinates) of docking 800 naturally occurring small molecules to the active site of a human oxidoreductase enzyme, the type 1 17-beta hydroxysteroid dehydrogenase (structure coordinates from PDB's 1equ).

Assume that we are only interested in the subset of steroid-like molecules from the original set we docked to 1equ. In other words, we want to extract from the original file only those molecules that have a certain substructure pattern (the steroid scaffold). Babel allows us to use SMARTS patterns to specify a substructure that we want either present or absent from the final set. Find out how to do this in babel, and use a SMARTS pattern that represents the basic steroid skeleton.
+ Hints: find the skeleton (the 4 rings only) in the top right diagram in wikipedia ([Wikipedia - Steroid](http://en.wikipedia.org/wiki/Steroid)), draw it on the editor available online from Xemistry (for simplicity make sure you use only the four fused rings), and copy the SMILES string from the box above the editor. Remember that SMILES strings are valid SMARTS expressions! Now you have a SMARTS pattern you can use to get a subset of the steroids using babel. How many molecules did you get? Save this new file as _1equ_ligonly_smarts.sdf_.

Visualise (with __Chimera__) your results to check that they indeed contain the SMARTS pattern you searched for. Chimera does not appear to be splitting the structures into models so you may want to try out the command
```
--splitinto
```
in babel to get the results in separate files. You can also use the -m command (which is supposed to print out multiple files), but it seems to be working only after the SMARTS subset has been created (at least for me!).

## Ready-to-use openbabel programs - Exploring other functions
There are many other things OpenBabel can do, although what is available through the ready-made programs is definitely restricted. Here we will explore a couple more useful functions.

Calculate some basic properties of the molecules in the original dataset (1equ_ligonly_edited.sdf) using the __obprop__ program. You will probably get a series of warnings. Without accessing the actual source code, it is difficult to know exactly what the problem is, but this is not unusual. As there are some online property calculators on the web, you could use the SMILES strings produced with obprop to submit these molecules online and check what results other programs are giving for a randomly selected molecule from your set (e.g. use the PubChem Identity search to retrieve the molecule of interest and then check its properties, or use the MolInspiration property calculator and paste in your SMILES string).

Use __obfit__ to superimpose all your steroid-containing molecules on the first entry of your file containing them (or to make it easier to see the result, you can just isolate two molecules and superimpose one on the other). Use any SMARTS pattern you want (suggestion: try the five-membered ring of the steroid skeleton). If you can, visualise the result (with Chimera, pymol or similar program). You can imagine how such a superposition might be useful for structure-activity studies. Obviously in this case the steroids were not in random orientations to begin with as they were the results of a docking simulation.

_End of practical_
