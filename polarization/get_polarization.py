#!/usr/bin/env python3
"""
Takes a list of directories containing VASP OUTCAR and POSCAR files, each
interpolating from a non-polar to a polar structure.
Prints out the polarization, by the modern theory of polarization.

Also supports a point-charge model, for 18-crown-6@Ba MnBr4.

Copyright 2025, Bernard Field
"""
import os
import argparse
import collections

import numpy as np
import matplotlib.pyplot as plt

from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.analysis.ferroelectricity import polarization


def get_polarity(dirs:list) -> polarization.Polarization:
    outcar_paths = [os.path.join(d, "OUTCAR") for d in dirs]
    poscar_paths = [os.path.join(d, "POSCAR") for d in dirs]
    # Convert OUTCAR to OUTCAR.gz as required
    outcar_paths = [p if os.path.isfile(p) else p+'.gz' for p in outcar_paths]

    structures = [Structure.from_file(p) for p in poscar_paths]
    outcars = [Outcar(o) for o in outcar_paths]

    polarity = polarization.Polarization.from_outcars_and_structures(outcars = outcars,
                                                                     structures = structures, 
                                                                     calc_ionic_from_zval = True)
    return polarity

def get_formal_polarity(dirs:list, zval_dict=dict(Br=-1,Mn=2,Ba=2,O=0,C=0,H=0)) -> polarization.Polarization:

    poscar_paths = [os.path.join(d, "POSCAR") for d in dirs]
    structures = [Structure.from_file(p) for p in poscar_paths]
    p_ions = [polarization.get_total_ionic_dipole(stru, zval_dict) for stru in structures]
    p_elecs = [[0.,0.,0.]] * len(p_ions)
    return polarization.Polarization(p_elecs, p_ions, structures)


def get_data_from_polarity(polarity) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Outputs:
        same-branch polarization (muC per cm2)
        Polarization change (inb basis of polar lattice vectors)
        Polarization change (in microC/cm^2)
    """
    same_branch = polarity.get_same_branch_polarization_data(convert_to_muC_per_cm2=True, all_in_polar=True)
    pol_change = polarity.get_polarization_change()
    pol_change_norm = polarity.get_polarization_change_norm()
    return same_branch, pol_change, pol_change_norm


def plot(polarity):
    same_branch = polarity.get_same_branch_polarization_data(convert_to_muC_per_cm2=True, all_in_polar=True)
    # Lattice unit vectors
    A = polarity.structures[-1].lattice.matrix
    # Divide each unit vector by its length. (Reshape to broadcast along correct axis)
    A = A / np.linalg.norm(A, axis=1).reshape((3,1))
    # same_branch has components along the lattice vectors
    # (read Polarization.get_polarization_change_norm)
    y = np.linalg.norm(same_branch @ A, axis=1)
    plt.plot(y, 'o-')
    plt.ylabel(r'Polarization ($\mu$C/cm$^2$)')
    plt.xlabel('Step')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dirs', nargs='+', help="Directories containing polarization DFT data.")
    parser.add_argument('-i','--no-ionic', action='store_true',
                        help="Neglect the ionic component of polarization.")
    parser.add_argument('-e','--no-elec', action='store_true',
                        help="Neglect the electronic component of the polarization.")
    parser.add_argument('-f','--formal', action='store_true',
                        help="Ignore usual calculation. Instead use formal ionic charges.")
    parser.add_argument('-d','--decomposed', action='store_true',
                        help="Also get decomposed atomic formal polarities.")
    parser.add_argument('-p','--plot', action='store_true',
                        help="Plot evolution of polarisation.")
    args = parser.parse_args()
    
    if args.formal:
        polarity = get_formal_polarity(args.dirs)
    else:
        polarity = get_polarity(args.dirs)
    if args.no_ionic:
        polarity.p_ions = np.zeros(polarity.p_ions.shape)
    if args.no_elec:
        polarity.p_elecs = np.zeros(polarity.p_elecs.shape)

    same_branch, pol_change, pol_change_norm = get_data_from_polarity(polarity)

    print("Same-branch polarization")
    print(same_branch)

    print("Polarization change (in basis of polar lattice vectors)")
    print(pol_change)

    print("Magnitude of polarization change in microC/cm^2")
    print(pol_change_norm)
    
    if args.plot:
        plot(polarity)
    
    if args.decomposed:
        zval_dict_original = dict(Br=-1,Mn=2,Ba=2,O=0,C=0,H=0)
        for k,v in (('Br',-1), ('Mn', 2), ('Ba', 2)):
            zval_dict = collections.defaultdict(lambda: 0)
            zval_dict[k] = v
            polarity = get_formal_polarity(args.dirs, zval_dict)
            print(f"Partial polarization for {k}")
            same_branch, pol_change, pol_change_norm = get_data_from_polarity(polarity)
            print("Polarization change (in basis of polar lattice vectors)")
            print(pol_change)

            print("Magnitude of polarization change in microC/cm^2")
            print(pol_change_norm)

