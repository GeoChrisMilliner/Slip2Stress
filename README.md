# Slip2Stress
MATLAB scripts to invert fault slip data for 3D deviatoric stress tensor

README for stress inversion 

Slip2Stress

v.1.02

Reference: 
Milliner et al. (2022)
Fault Friction Derived from Fault Bend Influence on Coseismic Slip During the 2019 Ridgecrest Earthquake 
JGR


Description
Slip2Stress is a MATLAB package to invert fault slip data (with a given strike, dip and rake) to invert for a 3D stress tensor. 
The inversion method can be determinded by the user ranging from an iterative L1 norm, damped L2 norm or a conjugate gradient descent
The inversion is based on Wallce-Bott assumption and method outlined by Michael (1984, 1987).
Damping method is based upon Hardebeck and Michael (2006).
The package can be run as-is, using example slip data from the 2019 Ridgecrest earthquake sequence 


Copyright

The code can be freely used for research purposes only. 
In the case of publishing the results obtained by this code, 
please, refer to the paper of Milliner et al. (2022). If you intend 
to use the code for commercial purposes, you should contact 
the author (geomilliner@gmail.com) for providing with the commercial licence. 
The use of the software for commercial purposes with no commercial licence is prohibited.


Other references

Hardebeck, J.L. and Michael, A.J., 2006. Damped regional‐scale stress inversions: Methodology and examples for southern California and the Coalinga aftershock sequence. Journal of Geophysical Research: Solid Earth, 111(B11).
Michael, A.J., 1984. Determination of stress from slip data: Faults and folds, J. Geophys. Res. 89, 11.517-11.526.
Michael, A.J., 1987. Use of focal mechanisms to determine stress: A control study, J. Geophys. Res. 92, 357-368.
Vavryčuk, V., 2014. Iterative joint inversion for stress and fault orientations from focal mechanisms, Geophysical Journal International, 199, 69-77, doi: 10.1093/gji/ggu224.
