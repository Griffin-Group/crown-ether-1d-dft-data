#!/usr/bin/env python
"""
Symmetrises the MnBr4 tetrahedron, restoring full tetrahedral symmetry.

Core functions written with assistance from LBL IT Division's CBorg Coder
(modified GPT OSS 120b). Bernard Field designed the high-level algorithm,
debugged the code, and formatted it properly for a script.
"""

import argparse

import numpy as np

from pymatgen.core import Structure

def symmetrise_MnBr(stru:Structure):
    """
    Symmetrise the MnBr4 tetrahedron of my specific structure.
    Edits Structure in place.
    """
    # Grab the Mn, center
    Mn = list(filter(lambda x: str(x.specie) == 'Mn', stru))
    Mn = Mn[0].coords
    # Grab the Br coords
    Brs = list(filter(lambda x: str(x.specie) == 'Br', stru))
    Brs = np.array([x.coords for x in Brs])
    # Find the one pointing towards (-1-1-1), heuristically closest to origin.
    idx = np.argmin(np.linalg.norm(Brs, axis=1))
    # Symmetrise the tetrahedron
    new_Brs = align_tetrahedron(Brs, Mn, idx)
    # Set the new coordinates.
    i = 0
    for site in stru:
        if str(site.specie) == 'Br':
            site.coords = new_Brs[i]
            i += 1
    # Done

def rotation_matrix(axis, theta):
    """Create a rotation matrix for rotating around an axis by an angle theta."""
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([
        [aa + bb - cc - dd, 2*(bc + ad), 2*(bd - ac)],
        [2*(bc - ad), aa + cc - bb - dd, 2*(cd + ab)],
        [2*(bd + ac), 2*(cd - ab), aa + dd - bb - cc]
    ])

def align_tetrahedron(corners, center, designated_corner_idx=0):
    # Step 1: Define the ideal tetrahedron coordinates
    ideal_corners = np.array([
        [0, 0, 1],
        [np.sqrt(8/9), 0, -1/3],
        [-np.sqrt(2/9), np.sqrt(2/3), -1/3],
        [-np.sqrt(2/9), -np.sqrt(2/3), -1/3]
    ])

    # Step 2: Rotate to align with the designated corner
    designated_corner = corners[designated_corner_idx] - center
    z_axis = designated_corner / np.linalg.norm(designated_corner)

    # Find the rotation axis that will align the ideal [0,0,1] with z_axis
    ideal_z = np.array([0, 0, 1])
    rotation_axis = np.cross(ideal_z, z_axis)
    if np.allclose(rotation_axis, 0):  # Already aligned
        theta = 0
        rotation_axis = ideal_z
    else:
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        cos_theta = np.dot(ideal_z, z_axis)
        theta = np.arccos(cos_theta)

    R1 = rotation_matrix(rotation_axis, theta)
    rotated_ideal = np.dot(R1, ideal_corners.T).T   

    # Step 3: Rotate around the designated corner axis to minimize distances
    # Pick one of the non-designated corners
    non_designated_corners = corners[[i for i in range(len(corners)) if i != designated_corner_idx]] - center
    non_designated_ideal = rotated_ideal[1:]

    # Find the best rotation angle for the first non-designated corner
    target_corner = non_designated_corners[0]
    ideal_corner = non_designated_ideal[0]

    # Project out the designated corner direction
    target_projected = target_corner - np.dot(target_corner, z_axis) * z_axis
    ideal_projected = ideal_corner - np.dot(ideal_corner, z_axis) * z_axis

    # Calculate the angle between the projected vectors
    cos_angle = np.dot(target_projected, ideal_projected) / (np.linalg.norm(target_projected) * np.linalg.norm(ideal_projected))
    sin_angle = np.sqrt(1 - cos_angle**2)

    # Determine the rotation direction
    cross_product = np.cross(target_projected, ideal_projected)
    if np.dot(cross_product, z_axis) < 0:
        sin_angle = -sin_angle

    # Calculate the rotation angle
    rotation_angle = np.arctan2(sin_angle, cos_angle)

    # Apply the opposite rotation to undo it.
    R2 = rotation_matrix(z_axis, -rotation_angle)
    final_rotated_ideal = np.dot(R2, rotated_ideal.T).T

    # Check if swapping the other two ideal corners improves alignment
    remaining_ideal = final_rotated_ideal[2:]
    remaining_targets = non_designated_corners[1:]

    # Calculate distances before swapping
    dist_before = np.sum((remaining_ideal - remaining_targets)**2)

    # Try swapping
    swapped_ideal = np.array([remaining_ideal[1], remaining_ideal[0]])
    dist_after = np.sum((swapped_ideal - remaining_targets)**2)

    # Use the swap if it improves alignment
    if dist_after < dist_before:
        final_rotated_ideal[2:] = swapped_ideal

    # Step 4: Scale to match the average bond length
    avg_bond_length = np.mean(np.linalg.norm(corners - center, axis=1))
    scaling_factor = avg_bond_length
    scaled_ideal = final_rotated_ideal * scaling_factor

    # Step 5: Translate to the center position
    transformed_ideal = scaled_ideal + center

    return transformed_ideal

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    args = parser.parse_args()

    stru = Structure.from_file(args.input)
    symmetrise_MnBr(stru)
    stru.to_file(args.output)
