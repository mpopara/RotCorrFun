# RotCorrFun
Python script to compute rotational correlation function of a protein backbone N-H vector 

## General description

NMR and fluorescence spectroscopy are examples of experimental methods which are sensitive to inter- and intramolecular dynamics across range of timescales, from picoseconds to hours. However, only all-atom MD simulations can
provide atomistic insight into these motions. Scripts within this repository serve for calculation of rotational correlation function, _C_(&tau;<sub>c</sub>), of a vector of interest, using trajectories obtained from all-atom MD simulations. 
Reorientation of a vector is determined by both overall rotational diffusion of a protein and internal motions. Since global rotational diffusion is much slower compared to internal motions of backbone or side-chain vectors, 
rotational correlation function _C_(&tau;<sub>c</sub>) can be separated into correlation function for global and internal rotational motions<sup>1</sup>:

C(&tau;<sub>c</sub>) = C<sub>internal</sub>(&tau;<sub>c</sub>)_C_<sub>global</sub>(&tau;<sub>c</sub>) 

If MD trajectory is aligned, overall rotational diffusion of a protein will be removed, and only timescales of internal rotational motions can be determined from the analysis of correlation function.


Rotational correlation function can be expressed in terms of second Legendre polynomial of the dot product of the unit vector _**e**_ at time t and t+&tau;<sub>c</sub> 
The scalar product between two unit vectors equals to the cosine of the angle _&theta;_ between them, which further simplifies the equation:  

C(&tau;<sub>c</sub>) = &lang;P<sub>2</sub>(_**e**_(t)&sdot;_**e**_(t+&tau;<sub>c</sub>))&rang; = &lang;&frac12;(3cos<sup>2</sup>_&theta;_(t+&tau;<sub>c</sub>)-1)&rang;

Within this repository, two scripts are provided, namely _RotCorrFun_full_traj.py_ and _RotCorrFun_subtraj.py_. In the first case, rotational correlation function is calculated for the entire trajectory. In the latter case, 
trajectory is divided in _N_ subtrajectories of _subtraj_length_. The idea is that by subsequently averaging over several subtrajectories, less noisy _C_(&tau;<sub>c</sub>) can be obtained.
Length of subtrajectory should be fewfold longer than the overal rotational diffusion time of the protein (if its estimate is available from other source). 

Particularly convenient is to calculate rotational correlation function for the backbone N-H bond vector, because this allows to directly compare the results with NMR data of {<sup>1</sup>H}-<sup>15</sup>N labelled protein: 

```
vector_origin = "resname LYS and name N"
vector_end = "resname LYS and name H"
```
In this particular example, I studied reorientation of N-H bond vectors in lysine residues. Script calculates and saves as .txt file individual lysine correlation functions for each of the (sub)trajectories,
as well as the average correlation function across all lysine residues within a (sub)trajectory. Depending on the system under the study, if local environments for different residues vary, such averaging will likely not be meaningful.
Therefore, it is left to the user to decide if C(&tau;<sub>c</sub>) should be averaged over residues and/or trajectories.

Using _mdtraj_ text-based atom selection language, reorientation of any other vector of choice can be monitored by defining its origin and end, as given in the example above.

![Lys_N-H_rot_corr_fun](https://github.com/mpopara/RotCorrFun/assets/40856779/78424a1a-2f6b-406e-9b2e-c7fd37a3d758)

C(&tau;<sub>c</sub>) can be modelled as weighted sum of exponential decays, whose characteristic times correspond to rotational correlation times of internal and global rotational motions. 
Furthermore, amplitudes of exponential terms relate to so-called order parameter, S<sup>2</sup>, which is a measure of bond flexibility, and in the case of backbone vectors, it relates to secondary structure.
Full framework related to analysis of rotational correlation function using fluorescence, NMR and MD data is presented in work by Möckel _et al_, 2019.<sup>1</sup>

## Input file requirement

* time-ordered trajectory obtained from all-atom MD simulations, in any of the mdtraj compatible formats (dcd, xtc, nc,..). If global rotational diffusion of the protein is to be quantified, trajectory should not be aligned.
  Otherwise, if only internal motions of backbone and side-chain vectors are of interest, trajectory should be superimposed on the first frame in trajectory.
* topology as .pdb file
  

## Dependencies
_RotCorrFun_full_traj.py_ and _RotCorrFun_subtraj.py_ are Python scripts developed on Python 3.8.8. Scripts were tested under the following configuration:

* Windows 10
* Python 3.8.8
* numpy 1.23.0
* mdtraj 1.9.4


## References

1. Mockel, C.; Kubiak, J.; Schillinger, O.; Kuhnemuth, R.; Della Corte, D.; Schroder,
G. F.; Willbold, D.; Strodel, B.; Seidel, C. A. M.; Neudecker, P., Integrated NMR,
Fluorescence, and Molecular Dynamics Benchmark Study of Protein Mechanics and
Hydrodynamics. J Phys Chem B 2019, 123 (7), 1453-1480.

## Authors

* Milana Popara
