# System
import sys
import os
import functools

# Data
import pandas as pd

# Maths
import numpy as np

# Atomic toolkits
from ase.io import read, write
from rascal.representations import SphericalInvariants as SOAP
from rascal.neighbourlist.structure_manager import mask_center_atoms_by_species
import global_descriptor

# ML
from sklearn.svm import SVC, LinearSVC
from sklearn.model_selection import GridSearchCV
from sklearn.utils import shuffle

# Data import
def convert_to_ase_with_mask(path):
    
    s = read(path, format='proteindatabank')

    if s.cell.volume <= 0:
        cell_size = 1000.0
        cell = np.ones(3) * cell_size
        s.cell = cell

    s.center()
    s.wrap(eps=1.0E-12)

    ca_indices = np.nonzero(s.arrays['atomtypes'] == 'CA')[0]
    s.arrays['numbers'][ca_indices] = 54
    mask_center_atoms_by_species(s, species_select=['Xe'])

    return s

# Load label data (4940 in total; some indices are not included due to no such structures)
df = pd.read_csv("3E2H.csv", index_col=0)

# Get shuffled indices (80-20)
n_data = len(df)
n_train = n_data*8//10
n_test = n_data - n_train
indices = shuffle(df, random_state=1).iloc[:n_data].index

# Load structures
structures = []
for index in indices:
    
    # import data
    path = '3E2H/3E2H_'+ str(index)  +'.pdb'    
    structures.append(convert_to_ase_with_mask(path))

# Define input argument for SOAP transformation
soap_parameters = dict(soap_type='PowerSpectrum', 
                       interaction_cutoff=5.0, 
                       max_radial=9,
                       max_angular=6,
                       gaussian_sigma_type='Constant', 
                       gaussian_sigma_constant=0.3,
                       cutoff_smooth_width=0.3,
                       normalize=False)

# Transform the structures by SOAP
soap = SOAP(**soap_parameters)
soaps = soap.transform(structures).get_features(soap)
n_centers = [np.count_nonzero(s.arrays['center_atoms_mask']) for s in structures]
split_idxs = np.cumsum(n_centers)[0:-1]
soaps = np.split(soaps, split_idxs)

# Define training set
X_train = soaps[:n_train]
X_test = soaps[n_train:]
Y_train = df["label"][indices[:n_train]].values
Y_test = df["label"][indices[n_train:]].values

## Compute the kernels (linear)
X_train_gaussian = []
for X in X_train:
    X_train_gaussian.append(X.mean(axis=0))

X_test_gaussian = []
for X in X_test:
    X_test_gaussian.append(X.mean(axis=0))

#X_train_gaussian = np.mean(X_train, axis=1)
#X_test_gaussian = np.mean(X_test, axis=1)

# Grid search with 5-fold CV
parameters = {'kernel':['linear', 'rbf'], 
              'C':[1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12], 
              'gamma':['scale']}
svc = SVC()
grid = GridSearchCV(svc, parameters, cv=5, iid=True)
grid_result = grid.fit(X_train_gaussian, Y_train)

means = grid_result.cv_results_['mean_test_score']
stds = grid_result.cv_results_['std_test_score']
params = grid_result.cv_results_['params']
for mean, stdev, param in zip(means, stds, params):
    print("%f (%f) with: %r" % (mean, stdev, param))
    

# Apply the best model
svc = SVC(C=grid_result.best_params_['C'],
          gamma=grid_result.best_params_['gamma'],
          kernel=grid_result.best_params_['kernel'])
svc.fit(X_train_gaussian, Y_train)

predicted_labels_gaussian_train = svc.predict(X_train_gaussian)
predicted_labels_gaussian_test = svc.predict(X_test_gaussian)

print('Train accuracy: %1.2f' % (np.sum(predicted_labels_gaussian_train == Y_train) / len(Y_train)))
print('Test accuracy: %1.2f' % (np.sum(predicted_labels_gaussian_test == Y_test) / len(Y_test)))
