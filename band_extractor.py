#!/usr/bin/env python3
"""
Extract selected band energies with k-point information from vaspout.h5
"""

import py4vasp as pv
import spglib
import numpy as np
import os
import sys
import argparse
from collections import defaultdict

def extract_selected_bands(vaspout_dir=".", band_indices=None, output_file="selected_bands.txt"):
    """
    Extract band energies for selected bands and write to file with k-point information.
    Full first Brillouin Zone is reconstructed starting from the Irreducible Brillouin Zone using symmetry operations.
    
    Parameters:
    -----------
    vaspout_dir : str
        Directory containing vaspout.h5 file
    band_indices : list or None
        List of band indices to extract (1-based indexing, e.g., [1, 2, 5])
        If None, extracts all bands
    output_file : str
        Name of output file
    """

    # Check if vaspout.h5 exists
    vaspout_file = os.path.join(vaspout_dir, "vaspout.h5")
    if not os.path.exists(vaspout_file):
        print(f"ERROR: {vaspout_file} not found!")
        sys.exit(1)
    
    print("=" * 70)
    print(f"Build band energy table from {vaspout_file}")
    print("=" * 70)
    print(f"Reading file: {vaspout_file}")
    print("=" * 70)
    print(f"Recovering symmetry operations from {vaspout_file}") 
    print("=" * 70)
    
    # Load the calculation
    try:
        calc = pv.Calculation.from_path(vaspout_dir)
        print("Loaded vaspout.h5")
    except Exception as e:
        print(f"Error loading vaspout.h5: {e}")
        sys.exit(1)
    
    # Get band structure data
    print("Extracting band structure data...")
    try:
        band_data = calc.band.to_dict()
        kpoint_data = calc.kpoint.to_dict()
        print("Successfully extracted band and k-point data")
    except Exception as e:
        print(f"Error extracting data: {e}")
        sys.exit(1)
    
    # Extract eigenvalues and k-point distances
    eigenvalues = band_data['bands']
    k_values = band_data['kpoint_distances']

    # Check orientation and ensure correct shape (nbands, nkpoints)
    if eigenvalues.shape[0] == len(k_values):
        eigenvalues = eigenvalues.T
        print(f"Transposed bands: now {eigenvalues.shape[0]} bands, {eigenvalues.shape[1]} k-points")

    nbands_total, nkpoints = eigenvalues.shape
    print(f"Total bands available: {nbands_total}")
    print(f"Total k-points: {nkpoints}")
    
    # Determine which bands to extract
    if band_indices is None:
        # Extract all bands
        bands_to_extract = list(range(1, nbands_total + 1))
        print(f"Extracting all {nbands_total} bands")
    else:
        # Validate band indices (convert to 0-based internally)
        bands_to_extract = []
        for idx in band_indices:
            if 1 <= idx <= nbands_total:
                bands_to_extract.append(idx)
            else:
                print(f"  Warning: Band {idx} is out of range (1-{nbands_total}), skipping")
        
        if not bands_to_extract:
            print("ERROR: No valid band indices selected!")
            sys.exit(1)
        
        print(f"Extracting {len(bands_to_extract)} selected bands: {bands_to_extract}")
    
    # Get k-point coordinates and weights
    kpoint_coords = kpoint_data['coordinates']
    kpoint_weights = kpoint_data['weights']
    mesh_dims = kpoint_data.get('mesh', None)
    if mesh_dims is None:
        print("  Note: No mesh information available in kpoint_data")
        mesh_dims = [30,30,30]

    # Get band occupations
    try:
        occupations = band_data['occupations']
        print("Successfully extracted occupations")
    except Exception as e:
        print(f"Error extracting occupations: {e}")
        sys.exit(1)

    occupation_type = []
    for bands in bands_to_extract:
        if occupations[0][bands - 1] == 1:
            occupation_type.append(1)
        elif occupations[0][bands - 1] == 0:
            occupation_type.append(0)
        else:
            occupation_type.append(occupations[0][bands - 1]) # metallic or partial occupation

    # Get lattice structure (POSCAR)
    try:
        structure = calc.structure[-1].to_dict()
        print("Successfully extracted lattice structure")
    except Exception as e:
        print(f"Error extracting lattice structure: {e}")
        sys.exit(1)
    lattice = structure['lattice_vectors']
    positions = structure['positions']
    elements = structure['elements']

    # Convert element symbols to atomic numbers (e.g., 'Si' -> 14)
    # For a simple mapping, we use a basic dictionary. For production, consider
    # using the `ase` library or a more complete mapping.
    atomic_numbers = []
    symbol_to_number = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
    'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
    'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
    'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
    'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
    'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
    'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
    'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86
    }

    for element in elements:
        if element in symbol_to_number:
            atomic_numbers.append(symbol_to_number[element])
        else:
            # Fallback: if element not in mapping, use a placeholder or raise error
            print(f"Warning: Atomic number for {element} not found, using 0.")
            atomic_numbers.append(0)
    
    # Assemble the cell tuple for spglib
    cell = (lattice, positions, atomic_numbers)

    # 4. Determine spglib version and handle accordingly
    symprec = 1e-5  # Symmetry tolerance
    try:
        dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
        
        # Check if we got a dictionary (old version) or an object (new version)
        if dataset is None:
            print("Symmetry analysis failed - dataset is None")
            
        # Detect return type
        is_dict = isinstance(dataset, dict)
        
        print(f"\n--- Symmetry Analysis Result ---")
        
        if is_dict:
            # Old spglib version (<2.5.0) - dataset is a dictionary
            print("  (Using spglib < 2.5.0 - dictionary format)")
            print(f"  Space Group Number: {dataset.get('number', 'N/A')}")
            print(f"  International Symbol: {dataset.get('international', 'N/A')}")
            print(f"  Hall Symbol: {dataset.get('hall', 'N/A')}")
            print(f"  Hall Number: {dataset.get('hall_number', 'N/A')}")
            print(f"  Choice: {dataset.get('choice', 'N/A')}")
            print(f"  Crystal System: {dataset.get('crystal_system', 'N/A')}")
            print(f"  Point Group: {dataset.get('pointgroup', 'N/A')}")
            
            # Get Bravais lattice from crystal system if available
            crystal_system = dataset.get('crystal_system', 'unknown')
            sg_number = dataset.get('number', 0)
            sg_symbol = dataset.get('international', 'N/A')
            crystal_system, bravais_lattice = get_crystal_info_from_sgnumber(sg_number, sg_symbol)

            print(f"  Crystal System: {crystal_system}")
            print(f"  Bravais Lattice: {bravais_lattice}")
            print(f"\n  {lattice[0][0]} {lattice[0][1]} {lattice[0][2]}")
            print(f"\n  {lattice[1][0]} {lattice[1][1]} {lattice[1][2]}")
            print(f"\n  {lattice[2][0]} {lattice[2][1]} {lattice[2][2]}\n")

            rotations = dataset['rotations']
            translations = dataset['translations']
            num_sym = len(rotations)
            print(f"Found {num_sym} symmetry operations from spglib")
        else:
            # New spglib version (>=2.5.0) - dataset is a SpglibDataset object
            # Note: 'crystal_system' and 'pointgroup' are NOT directly available as attributes
            print("  (Using spglib >= 2.5.0 - object format)")
            
            # Available attributes in SpglibDataset
            print(f"  Space Group Number: {dataset.number if hasattr(dataset, 'number') else 'N/A'}")
            print(f"  International Symbol: {dataset.international if hasattr(dataset, 'international') else 'N/A'}")
            print(f"  Hall Symbol: {dataset.hall if hasattr(dataset, 'hall') else 'N/A'}")
            print(f"  Hall Number: {dataset.hall_number if hasattr(dataset, 'hall_number') else 'N/A'}")
            print(f"  Choice: {dataset.choice if hasattr(dataset, 'choice') else 'N/A'}")

            # For crystal system and point group, we need to derive them from space group number
            sg_number = dataset.number if hasattr(dataset, 'number') else 0
            sg_symbol = dataset.international if hasattr(dataset, 'international') else 'N/A'
            crystal_system, bravais_lattice = get_crystal_info_from_sgnumber(sg_number, sg_symbol)
        
            print(f"  Crystal System: {crystal_system}")
            print(f"  Bravais Lattice: {bravais_lattice}")
            print(f"\n  {lattice[0][0]} {lattice[0][1]} {lattice[0][2]}")
            print(f"\n  {lattice[1][0]} {lattice[1][1]} {lattice[1][2]}")
            print(f"\n  {lattice[2][0]} {lattice[2][1]} {lattice[2][2]}\n")
                    
            # Also show Wyckoff letters if available
            if hasattr(dataset, 'wyckoffs'):
                print(f"  Wyckoff Letters: {dataset.wyckoffs[:5]}..." if len(dataset.wyckoffs) > 5 else f"  Wyckoff Letters: {dataset.wyckoffs}")

            rotations = dataset.rotations
            translations = dataset.translations
            num_sym = len(rotations)
            print(f"Found {num_sym} symmetry operations from spglib")

                
    except Exception as e:
        print(f"Error during symmetry analysis: {e}")
        print(f"Type of dataset: {type(dataset)}")
        
        # Alternative: use get_spacegroup() as fallback
        print("\n  Attempting fallback with get_spacegroup()...")
        try:
            sg_symbol = spglib.get_spacegroup(cell, symprec=symprec)
            print(f"  Space Group (from get_spacegroup): {sg_symbol}")
        except Exception as e2:
            print(f"  Fallback also failed: {e2}")

    full_kpoints, full_eigenvalues, kpoint_to_ir_map = reconstruct_full_brillouin_zone(kpoint_coords, kpoint_weights, lattice, rotations, mesh_dims, eigenvalues)
    print("bbb")
    # ========================================
    # Write to file
    # ========================================
    print(f"\nWriting data to {output_file}...")

    try:
        with open(output_file, 'w') as f:
            # Header
            f.write("#" + "=" * 68 + "\n")
            f.write("# BAND ENERGIES FROM VASP CALCULATION\n")
            f.write("#" + "=" * 68 + "\n")
            f.write(f"# Source: {vaspout_file}\n")
            f.write(f"# Number of bands total: {nbands_total}\n")
            f.write(f"# Number of k-points: {nkpoints}\n")
            f.write(f"# Number of bands extracted: {len(bands_to_extract)}\n")
            f.write(f"# Extracted bands: {', '.join(map(str, bands_to_extract))}\n")
            f.write(f"# Fermi energy set to 0 eV\n")
            f.write(f"{bravais_lattice}\n")
            f.write(f"{len(bands_to_extract)}   {nkpoints}\n")
            f.write(f"{lattice[0][0]} {lattice[0][1]} {lattice[0][2]}\n")
            f.write(f"{lattice[1][0]} {lattice[1][1]} {lattice[1][2]}\n")
            f.write(f"{lattice[2][0]} {lattice[2][1]} {lattice[2][2]}\n")
            f.write(f"# bands type (1 valence 0 conduction):\n")
            # Write band occupation types in header
            f.write(" ")
            for band_idx in bands_to_extract:
                f.write(f"{occupation_type[bands_to_extract.index(band_idx)]}" + " " * 3)
            f.write("\n")
            f.write("#" + "=" * 68 + "\n")
            f.write(f"HEADER END\n")
            f.write("#" + "=" * 68 + "\n")
            # Column headers
            f.write("# ")
            f.write(f"{'kx':>10} {'ky':>10} {'kz':>10} ")
            for band_idx in bands_to_extract:
                f.write(f"{f'band {band_idx:02d} E(eV)':>15} ")
            f.write("\n")
            
            f.write("# " + "-" * 68 + "\n")

            # Write data for each k-point
            for ik in range(len(full_kpoints)):                
                # k-point coordinates (fractional)
                f.write(f"{full_kpoints[ik][0]:10.6f} ")
                f.write(f"{full_kpoints[ik][1]:10.6f} ")
                f.write(f"{full_kpoints[ik][2]:10.6f} ")
                
                # Band energies for selected bands
                for band_idx in bands_to_extract:
                    # Convert to 0-based index
                    energy = full_eigenvalues[ik, band_idx - 1]
                    f.write(f"{energy:15.10f} ")
                f.write("\n")
            
            # Footer with statistics
            f.write("#" + "=" * 68 + "\n")
            f.write("# END OF FILE\n")
            f.write("#" + "=" * 68 + "\n")
        
        file_size = os.path.getsize(output_file)
        print(f"Successfully written to {output_file}")
        print(f"File size: {file_size:,} bytes")
        
    except Exception as e:
        print(f"  Error writing file: {e}")
        sys.exit(1)

    # Print summary
    print("\n" + "=" * 70)
    print("EXTRACTION SUMMARY")
    print("=" * 70)
    print(f"Extracted {len(bands_to_extract)} bands out of {nbands_total}")
    print(f"Extracted {nkpoints} k-points")
    print(f"Output file: {output_file}")
    print("\nEnergy range for extracted bands:")
    for band_idx in bands_to_extract:
        band_energies = eigenvalues[band_idx - 1, :]
        print(f" Band {band_idx:3d}: {band_energies.min():8.4f} eV to {band_energies.max():8.4f} eV")
    print("=" * 70)
    
    return eigenvalues[bands_to_extract[0]-1:bands_to_extract[-1], :], kpoint_coords

def get_crystal_info_from_sgnumber(sg_number, sg_symbol=None):
    """
    Derive crystal system, point group, and Bravais lattice from space group number.
    Based on the international tables of crystallography.
    """
    if 1 <= sg_number <= 2:
        return 'triclinic', 'Tr'
    elif 3 <= sg_number <= 15:
        if sg_symbol and 'P' in sg_symbol:
            return 'monoclinic', 'sM'
        elif sg_symbol and 'C' in sg_symbol:
            return 'monoclinic', 'oM'
        else:
            return 'monoclinic', 'sM'
    elif 16 <= sg_number <= 74:
        if(sg_symbol and 'P' in sg_symbol):
            return 'orthorhombic', 'sO'
        elif(sg_symbol and 'C' in sg_symbol):
            return 'orthorhombic', 'oO'
        elif(sg_symbol and 'I' in sg_symbol):
            return 'orthorhombic', 'bO'
        elif(sg_symbol and 'F' in sg_symbol):
            return 'orthorhombic', 'fO'
        else:
            return 'orthorhombic', 'sO'
    elif 75 <= sg_number <= 142:
        if sg_symbol and 'P' in sg_symbol:
            return 'tetragonal', 'sT'
        elif sg_symbol and 'I' in sg_symbol:
            return 'tetragonal', 'bT'
        else:
            return 'tetragonal', 'sT'
    elif 143 <= sg_number <= 167:
        if sg_symbol and 'P' in sg_symbol:
            return 'hexagonal', 'HEX'
        elif sg_symbol and 'R' in sg_symbol:
            return 'hexagonal', 'R'
        else:
            return 'hexagonal', 'HEX'
    elif 168 <= sg_number <= 194:
        return 'hexagonal', 'HEX'
    elif 195 <= sg_number <= 230:
        if sg_symbol and 'P' in sg_symbol:
            return 'cubic', 'sC'
        elif sg_symbol and 'I' in sg_symbol:
            return 'cubic', 'BCC'
        elif sg_symbol and 'F' in sg_symbol:
            return 'cubic', 'FCC'
        else:
            return 'cubic', 'sC'
    else:
        return 'unknown', 'unknown' #, 'unknown'

def reconstruct_full_brillouin_zone(ir_kpoints, ir_weights, lattice, rotations, mesh_dims, ir_eigenvalues, tolerance=1e-6):
    """
    Reconstruct the full Brillouin zone by applying symmetry operations to irreducible k-points.
    
    Parameters:
        ir_kpoints: list of irreducible k-points (fractional coordinates)
        ir_weights: normalized weights for irreducible k-points
        lattice: standard lattice for the given lattice type (real space)
        rotations: symmetry rotation matrices (from spglib)
        mesh_dims: tuple of (n1, n2, n3) for the full grid
        ir_eigenvalues: (m,n) matrix of energy eigenvalues with m number of irreducible k-points and n number of bands
        tolerance: tolerance for identifying duplicate k-points
    
    Returns:
        full_kpoints: list of all k-points in the full BZ (n_full, 3)
        full_energies: band energies for each full k-point (n_full, n_bands)
        kpoint_to_ir_map: mapping from full k-point index to irreducible index
    """

    ir_kpoints = np.array(ir_kpoints)
    ir_eigenvalues = np.array(ir_eigenvalues)
    ir_eigenvalues = np.transpose(ir_eigenvalues)
    rotations = np.array(rotations)
    reciprocal_lattice = np.linalg.inv(lattice)

    total_full_points = mesh_dims[0] * mesh_dims[1] * mesh_dims[2]
    integer_multiplicities = [int(round(w * total_full_points)) for w in ir_weights]

    full_kpoints = []
    kpoint_to_ir_map = []

    def fold_to_bz(k):
        """Fold k-point into the first Brillouin zone: range [-0.5, 0.5)."""
        return k - np.floor(k + 0.5)

    def is_duplicate(k, existing_kpoints):
        """Check if k is already in the list, accounting for BZ periodicity."""
        if len(existing_kpoints) == 0:
            return False, -1
        existing = np.array(existing_kpoints)
        diff = existing - k[np.newaxis, :]
        # Account for periodicity: shift by integers and check if difference is near zero
        #diff -= np.round(diff)
        dist = np.linalg.norm(diff, axis=1)
        idx = np.argmin(dist)
        if dist[idx] < tolerance:
            return True, idx
        return False, -1
    tolerance = 1e-6
    # Apply all symmetry rotations to each irreducible k-point
    for ir_idx, kpt in enumerate(ir_kpoints):
        # transform to cartesian coordinates
        k_basis = np.dot(reciprocal_lattice, kpt)
        num_point_diff = 0
        current_k_values = []
        for rot in rotations:
            # Rotate k-point: k' = R @ k (rotation in fractional coordinates)
            print(rot)
            k_rot = np.dot(rot, k_basis)
            duplicate, _ = is_duplicate(k_rot, current_k_values)
            if not duplicate:
                current_k_values.append(k_rot)
                k_rot = np.dot(lattice, k_rot)
                k_folded = fold_to_bz(k_rot)
                full_kpoints.append(k_folded)
                kpoint_to_ir_map.append(ir_idx)
                num_point_diff = num_point_diff + 1

        if(num_point_diff != integer_multiplicities[ir_idx]):
            print(f"Warning! Different number of k-points for {ir_idx} than expected. Expected {integer_multiplicities[ir_idx]}, found {num_point_diff}.")

    full_kpoints = np.array(full_kpoints)
    for i in range(len(full_kpoints)):
        full_kpoints[i] = np.dot(lattice, full_kpoints[i])
        full_kpoints[i] = fold_to_bz(full_kpoints[i])
    kpoint_to_ir_map = np.array(kpoint_to_ir_map)

    # Assign energies from irreducible k-point mapping
    full_energies = ir_eigenvalues[kpoint_to_ir_map]

    # Number of k-points in 1BZ vs number of k-points in k-mesh employed
    n_found = len(full_kpoints)
    print(f"Number of irreducible k-points: {len(ir_kpoints)}")
    print(f"Number of k-points used in the k-mesh: {total_full_points}")
    print(f"Number of k-points in full 1BZ: {n_found}")

    return full_kpoints, full_energies, kpoint_to_ir_map

def get_irreducible_kpoints_from_vaspout(vaspout_dir="."):
    """
    Extract irreducible k-points from vaspout.h5.
    """
    calc = pv.Calculation.from_path(vaspout_dir)
    kpoint_data = calc.kpoint.to_dict()
    
    ir_kpoints = kpoint_data['coordinates']  # Irreducible k-points
    ir_weights = kpoint_data['weights']      # Normalized weights
    mesh_dims = kpoint_data['mesh']          # Full mesh dimensions
    
    print(f"Number of irreducible k-points: {len(ir_kpoints)}")
    print(f"Sum of weights: {sum(ir_weights):.6f} (should be 1.0)")
    
    return ir_kpoints, ir_weights, mesh_dims

def find_little_group_operations(rotations, kpoint, tolerance=1e-6):
    """
    Find symmetry operations that leave a given k-point invariant (modulo reciprocal lattice).
    This is the "little group" of the k-point.
    """
    little_group = []
    
    for rot in rotations:
        # Transform k-point
        k_transformed = np.dot(rot, kpoint)
        # Check if k_transformed is equivalent to k (difference is integer)
        diff = k_transformed - kpoint
        diff = diff - np.round(diff)  # Bring to [-0.5, 0.5]
        
        if np.all(np.abs(diff) < tolerance):
            little_group.append(rot)
    
    return little_group

def save_symmetry_operations(rotations, translations, filename="symmetry_operations.txt"):
    """Save symmetry operations to a text file."""
    with open(filename, 'w') as f:
        f.write("# Space group symmetry operations\n")
        f.write(f"# Number of operations: {len(rotations)}\n")
        f.write("# Format: rotation_matrix (row-major) and translation vector\n\n")
        
        for i, (rot, trans) in enumerate(zip(rotations, translations)):
            f.write(f"Operation {i+1}:\n")
            f.write(f"  Rotation:\n")
            for row in rot:
                f.write(f"    [{row[0]:3d}, {row[1]:3d}, {row[2]:3d}]\n")
            f.write(f"  Translation: [{trans[0]:.6f}, {trans[1]:.6f}, {trans[2]:.6f}]\n\n")
    
    print(f"✓ Symmetry operations saved to {filename}")

def save_full_kpoints(full_kpoints, filename="full_kpoints.txt"):
    """Save reconstructed full Brillouin zone k-points to a file."""
    with open(filename, 'w') as f:
        f.write("# Full Brillouin zone k-points (fractional coordinates)\n")
        f.write(f"# Number of k-points: {len(full_kpoints)}\n")
        f.write("# kx, ky, kz\n")
        
        for kpt in full_kpoints:
            f.write(f"{kpt[0]:.6f}, {kpt[1]:.6f}, {kpt[2]:.6f}\n")
    
    print(f"✓ Full k-points saved to {filename}")

def main():
    """
    Command-line interface for band structure extraction.
    Main function to recover symmetries and 1BZ reconstruction.   
    """
    
    parser = argparse.ArgumentParser(
        description='Extract band structure from VASP vaspout.h5 file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Examples:
    # Extract bands 1-7 to a text file
    python extract_bands.py --bands 1-7 --output my_bands.txt
  
    # Extract specific bands 1,3,5,7 to CSV format
    python extract_bands.py --bands 1,3,5,7 --format csv --output bands.csv
  
    # Extract all bands with custom directory
    python extract_bands.py --dir /path/to/vasp/calculation --format txt
  
    # Extract bands near Fermi level (assuming around band 8)
    python extract_bands.py --bands 6-10 --output fermi_bands.txt
    """
    )
    
    # Input options
    parser.add_argument('--dir', '-d', type=str, default='.',
                        help='Directory containing vaspout.h5 (default: current directory)')
    
    parser.add_argument('--bands', '-b', type=str, default=None,
                        help='Band indices to extract. Examples: "1-7", "1,3,5,7", "1-5,8,10"')
    
    parser.add_argument('--output', '-o', type=str, default='selected_bands.txt',
                        help='Output file name (default: selected_bands.txt)')
    
    #parser.add_argument('--format', '-f', type=str, choices=['txt', 'csv'], default='txt',
    #                    help='Output format: txt (space-separated) or csv (comma-separated) (default: txt)')
    
    parser.add_argument('--all-bands', action='store_true',
                        help='Extract all bands (overrides --bands)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Parse band indices
    band_indices = None
    if args.all_bands:
        band_indices = None  # Will extract all bands
        print("Extracting all bands...")
    elif args.bands:
        band_indices = []
        parts = args.bands.split(',')
        for part in parts:
            part = part.strip()
            if '-' in part:
                # Handle range like "1-5"
                start, end = map(int, part.split('-'))
                band_indices.extend(range(start, end + 1))
            else:
                # Handle single number
                band_indices.append(int(part))
        print(f"Extracting bands: {band_indices}")
    else:
        # Default: bands 1-7
        band_indices = list(range(1, 8))
        print(f"No bands specified, using default: bands 1-7")
    
    # Text format (default)
    output_file = args.output
    if not output_file.endswith('.txt'):
        output_file = output_file.rsplit('.', 1)[0] + '.txt'
        
    extract_selected_bands(
        vaspout_dir=args.dir,
        band_indices=band_indices,
        output_file=output_file
    )

if __name__ == "__main__":
    main()