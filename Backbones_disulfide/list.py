import os

# Get all .pdb files in the current directory
pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]

# Write to a file, one name per line
with open('pdb_list.txt', 'w') as out_file:
    for pdb in pdb_files:
        out_file.write(f"{pdb}\n")

