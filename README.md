## Boltzmann distribution project
**Authors :** Maxime Borry, Ilyes Abdelhamid

### Requirements :  
- GROMACS (version > 5.0.7)
- Python (version > 2.7)
- R (version > 3)

### Usage :
`python boltzmann_project.py`

### Output :  
#### Files (in ./data/):
- `mass.txt` : mass of each atom of the pdb file, in g/mol
- `speed.csv` : speed components (x,y,z)of each atom, at each frame, 18699 atom per frame.
- `speed_vector.txt` : speed of each atom, at each frame, 18699 atom per frame.

#### Plots (in ./plot/) :
- `kinetik.png` : kinetik energy plot in function of time
- `speed_carbon_theo.png` : Boltzmann distribution of 1 speed component of the carbon, at 300K
- `temp_time.png` : temperature plot in function of time (calculated with equipartition and Maxwell-Boltzmann)
- prob_speed : PDF of speed x component for the carbon, on the first frame, with Maxwell-Boltzmann distribution, and fitted Gaussian function.

#### Other (in ./plot/) :
- `boltzmann.r.Rout` : R Output
- `boltzmann.RData` : R objects
