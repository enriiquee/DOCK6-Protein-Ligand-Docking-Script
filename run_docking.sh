#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Uso: $0 NOMBRE_PROTEINA"
    exit 1
fi

nombre_proteina="$1"
#############################################
#1.  Crear archivo INSPH:
#############################################

find . -name "temp*" -delete
rm -f rec.sph INSPH OUTSPH grid.in grid.nrg grid.cnt grid.bmp grid_tmp.in gridinfo.out *.sph

cd 002_surface_spheres
cat <<EOF > INSPH
../001_structure/${nombre_proteina}_rec_surface.dms
R
X
0.0
4.0
1.4
${nombre_proteina}_receptor_noH.sph
EOF

echo "DONE 1: Crear archivo INSPH"
# Eliminamos los archivos generados de otras ejecuciones


# Ejecutamos el comando necesasrio

sphgen -i INSPH -o OUTSPH

#############################################
#2. Seleccionar Esferas Relevantes
#############################################

sphere_selector ${nombre_proteina}_receptor_noH.sph ../001_structure/${nombre_proteina}_ligand_wH.mol2 10.0

echo "DONE 2. Seleccionar Esferas Relevantes"

#############################################
#3. Generar Caja 
#############################################

#Generamos archivo de configuración
cd ../003_gridbox

cat <<EOF > showbox.in
Y
8.0
../002_surface_spheres/selected_spheres.sph
1
${nombre_proteina}.box.pdb
EOF

#Arrancamos el programa
showbox < showbox.in

echo "DONE 3. Generar Caja"

#############################################
#4. Grid Generation
#############################################
#Generamos archivo de configuración Grid
cat <<EOF > grid.in
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
receptor_file                             ../001_structure/${nombre_proteina}_receptor_dockprep.mol2
box_file                                  ${nombre_proteina}.box.pdb
vdw_definition_file                       /home/eperez/dock6/parameters/vdw_AMBER_parm99.defn
score_grid_prefix                         grid
EOF


# Ejecuta el comando "grid" con el archivo de configuración actualizado
grid -i grid.in -o gridinfo.out

echo "DONE 4. Grid Generation"

#############################################
#5. Minimización Energia
#############################################

cd ../004_dock
cat <<EOF > min.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ../001_structure/${nombre_proteina}_ligand_wH.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes      
rmsd_reference_filename                                      ../001_structure/${nombre_proteina}_ligand_wH.mol2
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
ligand_outfile_prefix                                        ${nombre_proteina}.ligand.min
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

#Run dock
dock6 -i min.in -o min.out

echo "DONE 5. Minimización Energia"

########################
#Docking
########################

cd ../004_dock

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
ligand_atom_file                                             ./${nombre_proteina}.ligand.min_scored.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes  
rmsd_reference_filename                                      ./${nombre_proteina}.ligand.min_scored.mol2
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


echo "DONE: Docking" 

#Footprinting

cd 004_dock
cat <<EOF > footprint.in
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             ${nombre_proteina}.ligand.min_scored.mol2
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
fps_score_footprint_reference_mol2_filename                  ../001_structure/${nombre_proteina}_ligand_wH.mol2
fps_score_foot_compare_type                                  Euclidean
fps_score_normalize_foot                                     no
fps_score_foot_comp_all_residue                              yes
fps_score_receptor_filename                                  ../001_structure/${nombre_proteina}_receptor_dockprep.mol2
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