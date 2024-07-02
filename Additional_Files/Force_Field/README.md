NOTE: For convenience when packing the reverse micelles, the sodium is considered a part of the AOT molecule. It can be considered a separate molecule if you simply delete the relevant lines from each file. There are no bonds connecting the sodium to anything, so there is no need to search through the extra parameters in the ITP files to delete further entries. 

* "OPLS_AOT_Na.pdb" - generic structure file with atoms in the correct order for associated ITP files
* "OPLS-Std.itp" - the "standard" OPLS force field for AOT, using literature values to fill in missing parameters (see paper for details)
* "OPLS-CM5.itp" - OPLS-Std force field with updated partial charges on all atoms, using the CM5 partial charge calculation scheme (see paper for details)
* "OPLS-RESP.itp" - OPLS-Std force field with updated partial charges on all atoms, using the RESP partial charge calculation scheme (see paper for details)
