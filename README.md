# How-to-think-about-Shape

[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

Publicly available code and files supporting our work. Please cite the appropriate paper if you make use of these shape metrics.  

Development of shape metrics:

Gale, C.D.; Derakhshani-Molayousefi, M.; Levinger, N.E.; How to Characterize Amorphous Shapes: The Tale of a Reverse Micelle. *J Phys Chem B*, **2022**, *126* (4), 953-963. DOI: [10.1021/acs.jpcb.1c09439](https://doi.org/10.1021/acs.jpcb.1c09439)

Analysis of the shape of AOT reverse micelles and an OPLS parameterization of AOT:

Gale, C.D.; Derakhshani-Molayousefi, M.; Levinger, N.E.; Shape of AOT Reverse Micelles: The Mesoscopic Assembly is More Than the Sum of the Parts. *J Phys Chem B*, **2024**. DOI: [10.1021/acs.jpcb.4c02569](https://doi.org/10.1021/acs.jpcb.4c02569)

## Shape Analysis Code
We provide the code as is. To make use of the code, you will need to copy and paste the relevant functions into your own work. Thus, you may need to modify the code to suit your specific needs. For the most part, we expect any modifications will primarily be linked to the source of the data you are evaluating. We are chemists and this code is primarily designed for working with molecular dynamics (MD) simulations; any application outside of MD will likely require modification. For chemists performing MD, we used the GROMACS package and so had certain expectations of the file types. Make sure you update the code for your MD engine.

## Additional Files

We also provide several files which may be useful to people wishing to replicate or extend our simulations of AOT reverse micelles. These include:

* GROMACS style molecule ITP (and associated PDB) files for our OPLS parameterizations of AOT, including what we have called the OPLS-Std, OPLS-CM5, and OPLS-RESP force fields.
* 20 equilibrated starting structures from each simulation representing as many different unique shapes as possible, in PDB format. 

The diverse equilibrated structures were pulled from our 1 microsecond production run simulations, with stoichiometry of 42 AOT and Na+ ions, 210 water molecules, and 1500 isooctane solvent molecules (see 2024 paper for more details). We have generated a diverse set of configurations using K-means clustering based on the convexity and CPE values. While the data does not really fall into unique clusters, this method effectively partitions the shape space into N partitions of roughly equal size, which guarantees that the frames are as different from each other as possible. Note that running the same algorithm again will result in completely different sets of frames, but each calculation will still result in a diverse set of frames making the frames presented diverse although not unique. 

