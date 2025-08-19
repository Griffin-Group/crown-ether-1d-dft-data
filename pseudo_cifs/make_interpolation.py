#!/usr/bin/env python3

from pymatgen.core import Structure

if __name__ == "__main__":
    polar = Structure.from_file('symmetrised_original_primitive.vasp')
    nonpolar = Structure.from_file('pseudo_primitive_square_sorted.vasp')
    interp_list = nonpolar.interpolate(polar, interpolate_lattices=True, nimages=10)
    for i in range(len(interp_list)):
        interp_list[i].to_file(f'interpolation/POSCAR_{i}.vasp')
