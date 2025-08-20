# Data for DFT polarization calculations in "Supramolecular assembly of metal halide molecular wires"

If you use this data, please cite:

	Heqing Zhu, Cheng Zhu, Han K.D. Le, Daniel Cabeda, Bernard Field, Chuxi Wen, Alexander M. Oddo, Yuxin Jiang, Lihini Jayasinghe, Yu Shan, Lior Verbitsky, Harishankar Jayakumar, Sinéad M. Griffin, Eran Rabani, & Peidong Yang, "Supramolecular assembly of metal halide molecular wires", Nature Chemistry 2025 (under review).

This dataset is the data for the density functional theory (DFT) calculations used to calculate the polarization of (18-crown-6@Ba)MnBr4.
It was created by Bernard Field (Lawrence Berkeley National Laboratory).
Questions related to this data specifically should be directed to Bernard Field or Sinéad Griffin (Lawrence Berkeley National Laboratory).

## Computational info

DFT calculations were performed in VASP version 6.4.3 using standard settings on NERSC's Perlmutter supercomputer.

Polarization was calculated according to the Modern Theory of Polarization,
as implemented in VASP and Pymatgen.

The initial structure was provided by Heqing Zhu and should also be available on https://www.ccdc.cam.ac.uk/structures/ (see the paper associated with this work for reference numbers).

## Workflow

An initial structural relaxation (without magnetism) was performed in `1.relax`
followed by `2.relax2`
(where `2.relax2` has KPOINTS used in the single-point calculations).
The charge density (without spin) was calculated in `3.chg`, followed by the
magnetisation density in `4.mag`, each using the prior calculation as the
initial conditions.
We assumed a ferromagnetic configuration on Mn.
The final structure was obtained in `5.magrelax`, where a relaxation including
magnetism was performed.
A single-point calculation (with magnetism) was performed on this final structure in `6.mag2`.

The single-point calculations included
[Bader charge analysis](https://theory.cm.utexas.edu/henkelman/research/bader/),
although while somewhat helpful for diagnosing the oxidation states it was not
needed in the final paper.

Using the final structure, a non-polar reference structure was made in 
`pseudo_cifs/` (specifically, `pseudo_primitive_square_sorted.vasp`),
as well as the polar structures with symmetrised MnBr4 tetrahedra.
The polarization was then calculation in `polarization/`, where for each step an
initial single-point calculation for pre-converging the WAVECAR (`1.mag`) was
followed up by a self-consistent calculation of the polarization (`2.pol`)
(`LCALCPOL=.TRUE.`).

The polarization was then calculated using Pymatgen.

## Licenses

The contents are Copyright 2025 Bernard Field (as far as possible).

The data and documentation is licensed under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)
(CC-BY 4.0).

The code files (.py, slurm\_script, .sbatch) are licensed under the 
[MIT license](https://opensource.org/license/mit).

## Acknowledgements

We thank Ella Banyas for assistance and advice regarding the polarization calculations.

This research used resources of the National Energy Research Scientific
Computing Center, which is supported by the Office of Science of the U.S.
Department of Energy under Contract No. DE-AC02-05CH11231.
