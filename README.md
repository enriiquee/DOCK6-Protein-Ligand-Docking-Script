# DOCK6 Protein-Ligand Docking Script

This script streamlines the protein-ligand docking process using DOCK6.

## Prerequisites

1. Ensure you have DOCK6 installed and set up on your machine.
2. Modify the `base_directory` variable in the script to point to your desired base directory.
3. Ensure you have all required input files in the respective directories.

## Usage

To run the script:

```bash
./your_script_name.sh PROTEIN_NAME
```


Replace `your_script_name.sh` with the actual name of your script and `PROTEIN_NAME` with the name of your protein.

## Steps

1. **INSPH File Creation**: The script begins by generating an INSPH file.
2. **Select Relevant Spheres**: Using a given ligand file, the script selects spheres relevant for docking.
3. **Generate Box**: A box is generated for the docking process.
4. **Grid Generation**: Grids are generated for the docking calculations.
5. **Energy Minimization**: The ligand undergoes energy minimization.
6. **Docking**: The ligand is docked to the protein.
7. **Footprinting**: The docking footprints are saved.

After all the steps, the script cleans up the temporary files created during the process.

## Troubleshooting

If you encounter any issues, ensure that:
- All input files are correctly placed.
- DOCK6 parameters are correctly set.
- Necessary permissions are granted to the script.

## Contribution

Feel free to fork this repository, make changes, and submit pull requests. Any contributions, feedback, or issues are welcome!
