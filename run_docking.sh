#!/bin/bash

# Name: Enrique
# Email: enrique.perez@ucm.es

# Check protein name. 
if [ "$#" -ne 1 ] || [ -z "$1" ]; then
    echo "Usage: $0 PROTEIN_NAME" >&2
    exit 1
fi

protein_name="$1"
base_directory="$(pwd)"
dock6_parameters="/home/eperez/dock6/parameters"

# Create necessary directories if they don't exist
mkdir -p "$base_directory"/{002_surface_spheres,003_gridbox,004_dock,006_virtual_screen_mpi,007_cartesian_min,008_rescore}

# Function to clean up files
cleanup() {
    find . -name "temp*" -delete
    rm -f rec.sph INSPH OUTSPH grid.in grid.nrg grid.cnt grid.bmp grid_tmp.in gridinfo.out *.sph
}



#################################################
# 1. Create INSPH File
#################################################
cd "$base_directory/002_surface_spheres" || exit
cleanup
cat > INSPH <<EOF
../001_structure/${protein_name}_rec_surface.dms
R
X
0.0
4.0
1.4
${protein_name}_receptor_noH.sph
EOF

echo "DONE 1: Crear archivo INSPH"

# Running the necessary command

sphgen -i INSPH -o OUTSPH

#############################################
#2. Seleccionar Esferas Relevantes
#############################################

sphere_selector ${protein_name}_receptor_noH.sph ../001_structure/${protein_name}_ligand_wH.mol2 6.0
echo "DONE 2: Select Relevant Spheres"

#############################################
# 3. Generate Box
#############################################

#Generamos archivo de configuración
cd "$base_directory/003_gridbox" || exit

cat > showbox.in <<EOF
Y
8.0
../002_surface_spheres/selected_spheres.sph
1
${protein_name}.box.pdb
EOF

#Arrancamos el programa
showbox < showbox.in

echo "DONE 3: Generated Box"

#############################################
#4. Grid Generation
#############################################
#Generamos archivo de configuración Grid
cat > grid.in <<EOF
compute_grids                             yes
grid_spacing                              0.4
output_molecule                           no
contact_score                             no
energy_score                              yes
energy_cutoff_distance                    9999
atom_model                                a
attractive_exponent                       6
repulsive_exponent                        9
distance_dielectric                       yes
dielectric_factor                         4
allow_non_integral_charges                yes
bump_filter                               yes
bump_overlap                              0.75
receptor_file                             ../001_structure/${protein_name}_receptor_dockprep.mol2
box_file                                  ${protein_name}.box.pdb
vdw_definition_file                       /home/eperez/dock6/parameters/vdw_AMBER_parm99.defn
score_grid_prefix                         grid
EOF


# Execute the "grid" command with the configuration file updated
grid -i grid.in -o gridinfo.out

echo "============================================="
echo "Contenido de gridinfo.out"
echo "============================================="
cat gridinfo.out
echo "============================================="

echo "DONE 4. Grid Generation"

#############################################
# 5. Energy Minimization
#############################################

cd "$base_directory/004_dock" || exit
cat <<EOF > min.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ../001_structure/${protein_name}_ligand_wH.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes      
rmsd_reference_filename                                      ../001_structure/${protein_name}_ligand_wH.mol2
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ../003_gridbox/grid 
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_random_seed                                          0
simplex_restraint_min                                        yes
simplex_coefficient_restraint                                10.0
atom_model                                                   all
vdw_defn_file                                                /home/eperez/dock6/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /home/eperez/dock6/parameters/flex.defn
flex_drive_file                                              /home/eperez/dock6/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ${protein_name}.ligand.min
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

#Run dock
dock6 -i min.in -o min.out

echo "============================================="
echo "Contenido de min.out"
echo "============================================="

cat min.out

echo "============================================="


echo "DONE 5: Energy Minimization"

########################
# 6. Docking
########################

cd "$base_directory/004_dock" || exit
cat <<EOF > flex.in
conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               100.0
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            no
write_growth_tree                                            no
write_fragment_libraries                                     no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ./${protein_name}.ligand.min_scored.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes  
rmsd_reference_filename                                      ./${protein_name}.ligand.min_scored.mol2
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ../002_surface_spheres/selected_spheres.sph
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ../003_gridbox/grid
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                /home/eperez/dock6/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /home/eperez/dock6/parameters/flex.defn
flex_drive_file                                              /home/eperez/dock6/parameters/flex_drive.tbl
ligand_outfile_prefix                                        flex.out
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

#Run docking
dock6 -i flex.in -o flex.out

echo "============================================="
echo "Contenido de flex.out"
echo "============================================="
cat flex.out
echo "============================================="

echo "DONE: 6. Docking" 

#################################################
# 7. Footprinting
#################################################

cd "$base_directory/004_dock" || exit
cat <<EOF > footprint.in
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             ${protein_name}.ligand.min_scored.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     no
continuous_score_secondary                                   no
footprint_similarity_score_primary                           yes
footprint_similarity_score_secondary                         no
fps_score_use_footprint_reference_mol2                       yes
fps_score_footprint_reference_mol2_filename                  ../001_structure/${protein_name}_ligand_wH.mol2
fps_score_foot_compare_type                                  Euclidean
fps_score_normalize_foot                                     no
fps_score_foot_comp_all_residue                              yes
fps_score_receptor_filename                                  ../001_structure/${protein_name}_receptor_dockprep.mol2
fps_score_vdw_att_exp                                        6
fps_score_vdw_rep_exp                                        9
fps_score_vdw_rep_rad_scale                                  1
fps_score_use_distance_dependent_dielectric                  yes
fps_score_dielectric                                         4.0
fps_score_vdw_fp_scale                                       1
fps_score_es_fp_scale                                        1
fps_score_hb_fp_scale                                        0
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                /home/eperez/dock6/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /home/eperez/dock6/parameters/flex.defn
flex_drive_file                                              /home/eperez/dock6/parameters/flex_drive.tbl
ligand_outfile_prefix                                        footprint.out
write_footprints                                             yes
write_hbonds                                                 yes
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

#Run Docking

dock6 -i footprint.in -o footprint.out

echo "============================================="
echo "Contenido de footprint.out"
echo "============================================="

cat footprint.out
echo "============================================="

echo "DONE 7. Footprinting"

#################################################
# 8. Virtual Screen
#################################################
cd "$base_directory/005_virtual_screen" || exit
cat > virtual.in <<EOF
conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               100.0
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            no
write_growth_tree                                            no
write_fragment_libraries                                     no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             VS_library_25K.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no 
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ../002_surface_spheres/selected_spheres.sph
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no 
bump_filter                                                  no 
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ../003_gridbox/grid
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no 
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
minimize_flexible_growth_ramp                                no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0 
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                /home/eperez/dock6/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /home/eperez/dock6/parameters/flex.defn
flex_drive_file                                              /home/eperez/dock6/parameters/flex_drive.tbl
ligand_outfile_prefix                                        virtual.out
write_orientations                                           no 
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

# Run dock virtual
dock6 -i virtual.in -O virtual.out



#Change directory:
cd "$base_directory/006_virtual_screen_mpi" || exit

dock6 -i ../005_virtual_screen

echo "============================================="
echo "Contenido de virtual.out"
echo "============================================="
cat virtual.out
echo "============================================="






echo "DONE 8. Virtual Screen"

#################################################
# 9. Cartesian Energy Minimization
#################################################

cd "$base_directory/007_cartesian_min" || exit

cat > cart_min.in <<EOF
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ../006_virtual_screen_mpi/virtual.out_scored.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     yes
continuous_score_secondary                                   no
cont_score_rec_filename                                      ../001_structure/${protein_name}_receptor_dockprep.mol2
cont_score_att_exp                                           6
cont_score_rep_exp                                           12
cont_score_rep_rad_scale                                     1
cont_score_use_dist_dep_dielectric                           yes
cont_score_dielectric                                        4.0
cont_score_vdw_scale                                         1
cont_score_es_scale                                          1
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                /home/eperez/dock6/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /home/eperez/dock6/parameters/flex.defn
flex_drive_file                                              /home/eperez/dock6/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ${protein_name}.virtualscreen_cart.min
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

# Run cartesian

dock6 -i cart_min.in -o cart_min.virtual.mpi.out

# Añadir líneas de separación
echo "============================================="
echo "Contenido de cart_min.virtual.mpi.out:"
echo "============================================="

cat cart_min.virtual.mpi.out

echo "============================================="

echo "DONE 9. Cartesian Energy Minimization"


##################################################
# 10. Rescoring
##################################################

cd "$base_directory/008_rescore" || exit

cat > rescore.in <<EOF
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ../007_cartesian_min/virtualscreen_cart.min_scored.mol2 
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     no
continuous_score_secondary                                   no
footprint_similarity_score_primary                           no
footprint_similarity_score_secondary                         no
pharmacophore_score_primary                                  no
pharmacophore_score_secondary                                no
descriptor_score_primary                                     yes
descriptor_score_secondary                                   no
descriptor_use_grid_score                                    no
descriptor_use_multigrid_score                               no
descriptor_use_continuous_score                              no
descriptor_use_footprint_similarity                          yes
descriptor_use_pharmacophore_score                           yes
descriptor_use_tanimoto                                      yes
descriptor_use_hungarian                                     yes
descriptor_use_volume_overlap                                yes
descriptor_fps_score_use_footprint_reference_mol2            yes
descriptor_fps_score_footprint_reference_mol2_filename       ../004_dock/${protein_name}.ligand.min_scored.mol2
descriptor_fps_score_foot_compare_type                       Euclidean
descriptor_fps_score_normalize_foot                          no
descriptor_fps_score_foot_comp_all_residue                   yes
descriptor_fps_score_receptor_filename                       ../001_structure/${protein_name}_receptor_dockprep.mol2
descriptor_fps_score_vdw_att_exp                             6
descriptor_fps_score_vdw_rep_exp                             12
descriptor_fps_score_vdw_rep_rad_scale                       1
descriptor_fps_score_use_distance_dependent_dielectric       yes
descriptor_fps_score_dielectric                              4.0
descriptor_fps_score_vdw_fp_scale                            1
descriptor_fps_score_es_fp_scale                             1
descriptor_fps_score_hb_fp_scale                             0
descriptor_fms_score_use_ref_mol2                            yes
descriptor_fms_score_ref_mol2_filename                       ../004_dock/${protein_name}.ligand.min_scored.mol2
descriptor_fms_score_write_reference_pharmacophore_mol2      no
descriptor_fms_score_write_reference_pharmacophore_txt       no
descriptor_fms_score_write_candidate_pharmacophore           no
descriptor_fms_score_write_matched_pharmacophore             no
descriptor_fms_score_compare_type                            overlap
descriptor_fms_score_full_match                              yes
descriptor_fms_score_match_rate_weight                       5.0
descriptor_fms_score_match_dist_cutoff                       1.0
descriptor_fms_score_match_proj_cutoff                       0.7071
descriptor_fms_score_max_score                               20
descriptor_fingerprint_ref_filename                          ../004_dock/${protein_name}.ligand.min_scored.mol2
descriptor_hms_score_ref_filename                            ../004_dock/${protein_name}.ligand.min_scored.mol2
descriptor_hms_score_matching_coeff                          -5
descriptor_hms_score_rmsd_coeff                              1
descriptor_volume_score_reference_mol2_filename              ../004_dock/${protein_name}.ligand.min_scored.mol2
descriptor_volume_score_overlap_compute_method               analytical
descriptor_weight_fps_score                                  1
descriptor_weight_pharmacophore_score                        1
descriptor_weight_fingerprint_tanimoto                       -1
descriptor_weight_hms_score                                  1
descriptor_weight_volume_overlap_score                       -1
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                /home/eperez/dock6/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /home/eperez/dock6/parameters/flex.defn
flex_drive_file                                              /home/eperez/dock6/parameters/flex_drive.tbl
chem_defn_file                                               /home/eperez/dock6/parameters/chem.defn
pharmacophore_defn_file                                      /home/eperez/dock6/parameters/ph4.defn
ligand_outfile_prefix                                        descriptor.out
write_footprints                                             yes
write_hbonds                                                 yes
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

# Run
dock6 -i rescore.in -o rescore.out


echo "============================================="
echo "Contenido de rescore.out"
echo "============================================="

cat rescore.out

echo "============================================="
echo "DONE 10. Rescoring"