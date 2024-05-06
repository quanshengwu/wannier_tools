#!/usr/bin/env python

# Copyright (C) 2011 Atsushi Togo
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

#
# A python tool based on Phononpy to generate phonon tight-binding Hamiltonian
# whose format is similar to wannier90_hr.dat produced by Wannier90.
# Author: Changming Yue, Institute of Physics, CAS.
# email: yuechangming8@gmail.com
# Usage: python phonon_hr.py -d --dim="3 3 3" -c POSCAR-uinicell -p band.conf
# where dim="3 3 3" is just for a supercell with 3*3*3 size. You can set your
# own size as you want.
#

import sys
import os
import numpy as np
from phonopy import Phonopy, PhonopyGruneisen, PhonopyQHA, __version__
import phonopy.file_IO as file_IO
from phonopy.cui.settings import Settings, PhonopyConfParser
from phonopy.cui.show_symmetry import check_symmetry
from phonopy.cui.phonopy_argparse import get_parser
from phonopy.units import (Wien2kToTHz, AbinitToTHz, PwscfToTHz, ElkToTHz,
                           SiestaToTHz, VaspToTHz, Hartree, Bohr)
from phonopy.structure.cells import print_cell, determinant
from phonopy.structure.atoms import atom_data, symbol_map
from phonopy.interface import create_FORCE_SETS, read_crystal_structure
from fcmat2hr import get_phonon_hr
from phonopy.structure.cells import get_reduced_bases

phonopy_version = __version__

# AA is created at http://www.network-science.de/ascii/.

def print_author():
    print("""
    A python tool based on Phononpy to generate phonon tight-binding Hamiltonian
    whose format is similar to wannier90_hr.dat produced by Wannier90.
    Author: Changming Yue, Department of Physics, Southern University of Science and Technology.
    email: yuechangming8@gmail.com
    """)
    

def print_phononpy():
    print("""        _
  _ __ | |__   ___  _ __   ___   _ __  _   _  ____ _____
 | '_ \| '_ \ / _ \| '_ \ / _ \ | '_ \| | | |  ||  ||  \\\\
 | |_) | | | | (_) | | | | (_) || |_) | |_| |  ||  ||__//
 | .__/|_| |_|\___/|_| |_|\___(_) .__/ \__, |  ||  ||  \\\\
 |_|                            |_|    |___/   ||  ||__// """)

    print_author()

def print_version(version):
    print(('%s' % version).rjust(44))
    print('')

def print_end():
    print("""                 _
   ___ _ __   __| |
  / _ \ '_ \ / _` |
 |  __/ | | | (_| |
  \___|_| |_|\__,_|
""")

def print_error():
    print("""  ___ _ __ _ __ ___  _ __
 / _ \ '__| '__/ _ \| '__|
|  __/ |  | | | (_) | |
 \___|_|  |_|  \___/|_|
""")

def print_attention(attention_text):
    print("*" * 67)
    print(attention_text)
    print("*" * 67)
    print('')

def print_error_message(message):
    print('')
    print(message)

def file_exists(filename, log_level):
    if os.path.exists(filename):
        return True
    else:
        error_text = "%s not found." % filename
        print_error_message(error_text)
        if log_level > 0:
            print_error()
        sys.exit(1)

def print_cells(phonon, unitcell_filename):
    print("Crsytal structure is read from \'%s\'." % unitcell_filename)
    supercell = phonon.get_supercell()
    unitcell = phonon.get_unitcell()
    primitive = phonon.get_primitive()
    p2p_map = primitive.get_primitive_to_primitive_map()
    mapping = np.array(
        [p2p_map[x] for x in primitive.get_supercell_to_primitive_map()],
        dtype='intc')
    s_indep_atoms = phonon.get_symmetry().get_independent_atoms()
    p_indep_atoms = mapping[s_indep_atoms]
    if unitcell.get_number_of_atoms() == primitive.get_number_of_atoms():
        print("-" * 32 + " unit cell " + "-" * 33)
        print_cell(primitive, stars=p_indep_atoms)
    else:
        u2s_map = supercell.get_unitcell_to_supercell_map()
        print("-" * 30 + " primitive cell " + "-" * 30)
        print_cell(primitive, stars=p_indep_atoms)
        print("-" * 32 + " unit cell " +  "-" * 33) # 32 + 11 + 33 = 76
        u2u_map = supercell.get_unitcell_to_unitcell_map()
        u_indep_atoms = [u2u_map[x] for x in s_indep_atoms]
        print_cell(unitcell, mapping=mapping[u2s_map], stars=u_indep_atoms)
    print("-" * 32 + " super cell " + "-" * 32)
    print_cell(supercell, mapping=mapping, stars=s_indep_atoms)
    print("-" * 76)

def print_settings(settings):
    run_mode = settings.get_run_mode()
    if run_mode == 'band':
        print("Band structure mode xixix")
    if run_mode == 'mesh':
        print("Mesh sampling mode")
    if run_mode == 'band_mesh':
        print("Band structure and mesh sampling mode")
    if run_mode == 'anime':
        print("Animation mode")
    if run_mode == 'modulation':
        print("Modulation mode")
    if run_mode == 'irreps':
        print("Ir-representation mode")
    if run_mode == 'qpoints':
        if settings.get_write_dynamical_matrices():
            print("QPOINTS mode (dynamical matrices written out)")
        else:
            print("QPOINTS mode")
    if (run_mode == 'band' or
        run_mode == 'mesh' or
        run_mode == 'qpoints') and settings.get_is_group_velocity():
        gv_delta_q = settings.get_group_velocity_delta_q()
        if gv_delta_q is not None:
            print("  With group velocity calculation (dq=%3.1e)" % gv_delta_q)
        else:
            print('')
    if run_mode == 'displacements':
        print("Creating displacements")
        if not settings.get_is_plusminus_displacement() == 'auto':
            if settings.get_is_plusminus_displacement():
                print("  Plus Minus displacement: full plus minus directions")
            else:
                print("  Plus Minus displacement: only one direction")
        if not settings.get_is_diagonal_displacement():
            print("  Diagonal displacement: off")

    print("Settings:")
    if settings.get_is_nac():
        print("  Non-analytical term correction: on")
    if settings.get_fc_spg_symmetry():
        print("  Enforce space group symmetry to force constants: on")
    if settings.get_fc_symmetry_iteration() > 0:
        print("  Force constants symmetrization: %d times" %
              settings.get_fc_symmetry_iteration())
    if settings.get_lapack_solver():
        print("  Use Lapack solver via Lapacke: on")
    if run_mode == 'mesh' or run_mode == 'band_mesh':
        print("  Sampling mesh: %s" % np.array(settings.get_mesh()[0]))

    if (np.diag(np.diag(settings.get_supercell_matrix())) \
            - settings.get_supercell_matrix()).any():
        print("  Supercell matrix:")
        for v in settings.get_supercell_matrix():
            print("    %s" % v)
    else:
        print("  Supercell: %s" % np.diag(settings.get_supercell_matrix()))
    if settings.get_primitive_matrix() is not None:
        print("  Primitive axis:")
        for v in settings.get_primitive_matrix():
            print("    %s" % v)

###################
# Command options #
###################
parser = get_parser()
(options, args) = parser.parse_args()
option_list = parser.option_list

# Set log level
log_level = 1
if options.verbose:
    log_level = 2
if options.quiet or options.is_check_symmetry:
    log_level = 0
if options.loglevel is not None:
    log_level=options.loglevel

#
# Phonopy interface mode
#
# Physical units: energy,  distance,  atomic mass, force
# vasp          : eV,      Angstrom,  AMU,         eV/Angstrom
# wien2k        : Ry,      au(=borh), AMU,         mRy/au
# abinit        : hartree, au,        AMU,         eV/Angstrom
# elk           : hartree, au,        AMU,         hartree/au
# pwscf         : Ry,      au,        AMU,         Ry/au
# siesta        : eV,      au,        AMU,         eV/Angstroem
#
if options.wien2k_mode:
    interface_mode = 'wien2k'
    from phonopy.interface.wien2k import write_supercells_with_displacements
    factor = Wien2kToTHz
    nac_factor = 2000.0
    distance_to_A = Bohr
elif options.abinit_mode:
    interface_mode = 'abinit'
    from phonopy.interface.abinit import write_supercells_with_displacements
    factor = AbinitToTHz
    nac_factor = Hartree / Bohr
    distance_to_A = Bohr
elif options.pwscf_mode:
    interface_mode = 'pwscf'
    from phonopy.interface.pwscf import write_supercells_with_displacements
    factor = PwscfToTHz
    nac_factor = 2.0
    distance_to_A = Bohr
elif options.elk_mode:
    interface_mode = 'elk'
    from phonopy.interface.elk import write_supercells_with_displacements
    factor = ElkToTHz
    nac_factor = 1.0
    distance_to_A = Bohr
elif options.siesta_mode:
    interface_mode = 'siesta'
    from phonopy.interface.siesta import write_supercells_with_displacements
    factor = SiestaToTHz
    nac_factor = Hartree / Bohr
    distance_to_A = Bohr
else:
    if options.vasp_mode:
        interface_mode = 'vasp'
        distance_to_A = 1.0
    else:
        interface_mode = None
        distance_to_A = 1.0
    from phonopy.interface.vasp import write_supercells_with_displacements
    from phonopy.interface.vasp import create_FORCE_CONSTANTS
    factor = VaspToTHz
    nac_factor = Hartree * Bohr

# Show title
if log_level > 0:
    print_phononpy()
    print_version(phonopy_version)

# Create FORCE_SETS for VASP (-f or --force_sets)
if options.force_sets_mode or options.force_sets_zero_mode:
    file_exists('disp.yaml', log_level)
    for filename in args:
        file_exists(filename, log_level)
    error_num = create_FORCE_SETS(interface_mode,
                                  args,
                                  options.symprec,
                                  is_wien2k_p1=options.is_wien2k_p1,
                                  log_level=log_level)
    if log_level > 0:
        print_end()
    sys.exit(error_num)

# Create FORCE_CONSTANTS (--fc or --force_constants)
if options.force_constants_mode:
    if len(args) > 0:
        file_exists(args[0], log_level)
        error_num = create_FORCE_CONSTANTS(args[0],
                                           options,
                                           log_level)
    else:
        print_error_message("Please specify vasprun.xml.")
        error_num = 1

    if log_level > 0:
        print_end()
    sys.exit(error_num)

# Parse the setting file
if len(args) > 0:
    if file_exists(args[0], log_level):
        phonopy_conf = PhonopyConfParser(filename=args[0],
                                         options=options,
                                         option_list=option_list)
        settings = phonopy_conf.get_settings()
else:
    phonopy_conf = PhonopyConfParser(options=options,
                                     option_list=option_list)
    settings = phonopy_conf.get_settings()

if options.is_graph_save:
    import matplotlib
    matplotlib.use('Agg')

##################
# Initialization #
##################
unitcell, optional_structure_file_information = read_crystal_structure(
    filename=settings.get_cell_filename(),
    interface_mode=interface_mode,
    chemical_symbols=settings.get_chemical_symbols(),
    yaml_mode=settings.get_yaml_mode())
unitcell_filename = optional_structure_file_information[0]

if unitcell is None:
    print_error_message("Crystal structure file of %s could not be found." %
                        unitcell_filename)
    if log_level > 0:
        print_error()
    sys.exit(1)

# Check unit cell
if np.linalg.det(unitcell.get_cell()) < 0.0:
    print_error_message("Determinant of the lattice vector matrix "
                        "has to be positive.")
    if log_level > 0:
        print_end()
    sys.exit(0)

# Set magnetic moments
magmoms = settings.get_magnetic_moments()
if magmoms is not None:
    if len(magmoms) == unitcell.get_number_of_atoms():
        unitcell.set_magnetic_moments(magmoms)
    else:
        error_text = "Invalid MAGMOM setting"
        print_error_message(error_text)
        if log_level > 0:
            print_end()
        sys.exit(1)

# Check crystal symmetry and exit (--symmetry)
if options.is_check_symmetry:
    check_symmetry(unitcell,
                   primitive_axis=settings.get_primitive_matrix(),
                   symprec=options.symprec,
                   distance_to_A=distance_to_A,
                   phonopy_version=phonopy_version)
    if log_level > 0:
        print_end()
    sys.exit(0)

# Phonon calculation mode: Band, mesh, qpoints, etc
run_mode = settings.get_run_mode()

# --factor (overwrite default factor for calculators)
if settings.get_frequency_conversion_factor() is not None:
    factor = settings.get_frequency_conversion_factor()

# --amplitude
if settings.get_displacement_distance() is None:
    if (interface_mode == 'wien2k' or
        interface_mode == 'abinit' or
        interface_mode == 'elk' or
        interface_mode == 'pwscf' or
        interface_mode == 'siesta'):
        displacement_distance = 0.02
    else: # default or vasp
        displacement_distance = 0.01
else:
    displacement_distance = settings.get_displacement_distance()

# Supercell matrix
if settings.get_supercell_matrix() is None:
    print_error_message("Supercell matrix (DIM or --dim) is not found.")
    if log_level > 0:
        print_end()
    sys.exit(1)
num_atom = unitcell.get_number_of_atoms()
num_satom = determinant(settings.get_supercell_matrix()) * num_atom
if settings.get_is_force_constants() == 'read':
    if file_exists("FORCE_CONSTANTS", log_level):
        fc = file_IO.parse_FORCE_CONSTANTS(filename="FORCE_CONSTANTS")
        fc_filename = "FORCE_CONSTANTS"

    if log_level > 0:
        print("Force constants are read from %s." % fc_filename)

    if fc.shape[0] != num_satom:
        error_text = ("Number of atoms in supercell is not consistent with "
                      "the matrix shape of\nforce constants read from ")
        if settings.get_is_hdf5():
            error_text += "force_constants.hdf5.\n"
        else:
            error_text += "FORCE_CONSTANTS.\n"
        error_text += ("Please carefully check DIM, FORCE_CONSTANTS, "
                       "and %s.") % unitcell_filename
        print_error_message(error_text)
        if log_level > 0:
            print_end()
        sys.exit(1)

else:
    if log_level > 1:
        print_end()
        sys.exit(1)
    else:
        file_exists("FORCE_SETS", log_level)

phonon = Phonopy(unitcell,
                 settings.get_supercell_matrix(),
                 primitive_matrix=settings.get_primitive_matrix(),
                 factor=factor,
                 is_auto_displacements=False,
                 dynamical_matrix_decimals=settings.get_dm_decimals(),
                 force_constants_decimals=settings.get_fc_decimals(),
                 symprec=options.symprec,
                 is_symmetry=settings.get_is_symmetry(),
                 use_lapack_solver=settings.get_lapack_solver(),
                 log_level=log_level)

supercell = phonon.get_supercell()
primitive = phonon.get_primitive()

# Set atomic masses of primitive cell
if settings.get_masses() is not None:
    phonon.set_masses(settings.get_masses())

# Print cells
if log_level > 1:
    print_cells(phonon, unitcell_filename)

# Set force constants
if settings.get_is_force_constants() == 'read':
    phonon.set_force_constants(fc)

# Impose cutoff radius on force constants
cutoff_radius = settings.get_cutoff_radius()
if cutoff_radius:
    phonon.set_force_constants_zero_with_radius(cutoff_radius)

# Enforce space group symmetry to force constants
if settings.get_fc_spg_symmetry():
    if log_level > 0:
        print('')
        print("Force constants are symmetrized by space group operations.")
        print("This may take some time...")
    phonon.symmetrize_force_constants_by_space_group()
    file_IO.write_FORCE_CONSTANTS(phonon.get_force_constants(),
                                  filename='FORCE_CONSTANTS_SPG')
    if log_level > 0:
        print("Symmetrized force constants are written into "
              "FORCE_CONSTANTS_SPG.")

# Imporse translational invariance and index permulation symmetry to
# force constants
if settings.get_fc_symmetry_iteration() > 0:
    phonon.symmetrize_force_constants(settings.get_fc_symmetry_iteration())

# Write FORCE_CONSTANTS
if settings.get_is_force_constants() == "write":
    if settings.get_is_hdf5():
        file_IO.write_force_constants_to_hdf5(phonon.get_force_constants())
        if log_level > 0:
            print("Force constants are written into force_constants.hdf5.")
    else:
        file_IO.write_FORCE_CONSTANTS(phonon.get_force_constants())
        if log_level > 0:
            print("Force constants are written into FORCE_CONSTANTS.")

# Show the rotational invariance condition (just show!)
if settings.get_is_rotational_invariance():
    phonon.get_rotational_condition_of_fc()

# Atomic species without mass case
symbols_with_no_mass = []
if primitive.get_masses() is None:
    for s in primitive.get_chemical_symbols():
        if (atom_data[symbol_map[s]][3] is None and
            s not in symbols_with_no_mass):
            symbols_with_no_mass.append(s)
            print_error_message(
                "Atomic mass of \'%s\' is not implemented in phonopy." % s)
            print_error_message("MASS tag can be used to set atomic masses.")

if len(symbols_with_no_mass) > 0:
    if log_level > 0:
        print_end()
    sys.exit(1)


phonon.set_dynamical_matrix()
dmat = phonon._dynamical_matrix

# rescale fcmat by THZ**2
fcmat = dmat._force_constants  * factor**2 # FORCE_CONSTANTS
smallest_vectors = dmat._smallest_vectors
#mass = dmat._mass
mass = dmat._pcell.get_masses()
multi = dmat._multiplicity

reduced_bases = get_reduced_bases(supercell.get_cell(), options.symprec)
positions = np.dot(supercell.get_positions(), np.linalg.inv(reduced_bases))
#for pos in positions: pos -= np.rint(pos)

relative_scale = np.dot(reduced_bases, np.linalg.inv(primitive.get_cell()))
super_pos=np.zeros((num_satom,3),dtype = np.float64)
for i in range(num_satom):
   super_pos[i] = np.dot(positions[i] , relative_scale)

p2s_map = dmat._p2s_map = primitive.get_primitive_to_supercell_map()
s2p_map = dmat._s2p_map = primitive.get_supercell_to_primitive_map()


num_satom = supercell.get_number_of_atoms()
num_patom = primitive.get_number_of_atoms()

# get hr
print "Generating phononpy tight binding Hamiltonian phonopy_TB.dat ..."
print "Force constants are rescaled by multiplying THZ**2 and by dividing multiplicity and sqrt(mass_iatom*mass_jatom) !"
get_phonon_hr(fcmat,smallest_vectors,mass,multi,super_pos,p2s_map,s2p_map,num_satom,num_patom)
print "phonopy_TB.dat generated! "
# END
if log_level > 0:
    print " to print end now .."
    print_end()


