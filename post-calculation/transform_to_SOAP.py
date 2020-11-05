import sys
import os
import numpy as np
from ase.io import read, write
from rascal.representations import SphericalInvariants as SOAP
from rascal.neighbourlist.structure_manager import mask_center_atoms_by_species

structure = read('6zpt.pdb', format='proteindatabank')

#reformat the cell if there is not any.
#wrap the structure into the cell box with tolerance eps
if structure.cell.volume <= 0:
    cell_size = 1000.0
    cell = np.ones(3) * cell_size
    structure.cell = cell
structure.center()
structure.wrap(eps=1.0E-12)

#mask H atoms
mask_center_atoms_by_species(structure, species_blacklist=['H'])

#zoom in backbone CA
#change the selected CA to Xe (54) and choose them
#ca_indices = np.nonzero(structure.arrays['atomtypes'] == 'CA')[0]
#structure.arrays['numbers'][ca_indices] = 54
#mask_center_atoms_by_species(structure, species_select=['Xe'])

soap_parameters = dict(soap_type='PowerSpectrum', # Basis of spherical harmonics and radial functions
                       interaction_cutoff=5.0, # Cutoff distance for the local atomic environments
                       max_radial=9, # Max number of radial basis functions
                       max_angular=6, # Max angular momentum number
                       gaussian_sigma_type='Constant', # Type of Gaussian width to use in building the atom density
                       gaussian_sigma_constant=0.3, # Gaussian width for building the atom density
                       cutoff_smooth_width=0.3, # Cutoff transition width
                       normalize=False) # Whether to normalize the SOAP vectors to have norm 1 -- for SVM, we probably want to use unnormalized SOAP vectors

#create transforming object soap
#compute SOAP vector soaps
soap = SOAP(**soap_parameters)
soaps = soap.transform(structure).get_features(soap)
