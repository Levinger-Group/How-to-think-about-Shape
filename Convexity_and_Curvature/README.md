# Convexity and Curvature

Convexity and both mean and Gaussian curvature are both computed from the same Willard-Chandler surface of your object. The creation of that Willard-Chandler surface can be somewhat computationally expensive, so it's typically best to compute both values at the same time. For this reason, we thought it would be best to include essentially 2 copies of the code for each. In the Jupyter Notebook, we'll go into detail about how it works, but won't put the whole thing together into useable functions. In the python files, we will combine those components into usable code that can be run to compute these values as well as detail how to read/write the results from/to file. There are python files for both a serial version and a parallel version available. 

Currently, PyTim is required to run these analyses and does not function properly on Windows operating systems. We have used it and it works great on Linux systems and we would assume it works well on MacOS systems, although we haven't tested that either. 

CG, 2021

Sega, M.; Hantal, G.; Fabian, B.; Jedlovszky, P.; Pytim: a Python Package for the Interfacial Analysis of Molecular Simulations. *J. Comp. Chem.*, **2018**, *39*, 2118-2125.

Willard, A.P.; Chandler, D.; Instantaneous Liquid Interfaces. *J. Phys. Chem. B*, **2010**, *114*, 1954-1958.
