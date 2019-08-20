# PRM
![alt text](https://zenodo.org/badge/203188792.svg)

PRM (Positional Recurrence Map), is a Matlab code to post-process ab initio molecular dynamics trajectories obtained with VASP, in order to obtain information on the structure and the dynamics of the system and, notably, to calculate neutron-weighted spectra. 
See prm.m header and functions headers for a detailed description of the code capabilities.

It was originally developed to calculate nuclear density maps to separate the static (structural distortions) and dynamic contibutions to the atomic motion (see Piovano, Perrichon, Boehm, Johnson and Paulus, Phys. Chem. Chem. Phys., 2016,18, 17398-17403   https://doi.org/10.1039/C5CP06464C).
It has then been extended to compute neutron-weighted spectra by power spectral density (PSD) and structural averages (see Mazzei, Perrichon, Mancini, Wahnström, Malavasi, Parker, Börjesson, Karlsson, J. Mater. Chem. A, 2019,7, 7360-7372   http://dx.doi.org/10.1039/c8ta06202a).
The latest version can also calculate the overtones and combination modes of the neutron spectra (see Perrichon, Jiménez-Ruiz, Mazzei, Rahman and Karlsson, J. Mater. Chem. A, 2019,7, 17626-17636, https://doi.org/10.1039/C9TA04056K).
