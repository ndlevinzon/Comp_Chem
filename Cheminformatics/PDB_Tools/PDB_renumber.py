# ND Levinzon, FAU Erlangen 2022
#  __     ________ _   _ _______
#  \ \   / /  ____| \ | |__   __|/\
#   \ \_/ /| |__  |  \| |  | |  /  \
#    \   / |  __| | . ` |  | | / /\ \
#     | |  | |____| |\  |  | |/ ____ \
#     |_|  |______|_| \_|  |_/_/    \_\


# Insert the name of the PDB file in between the ""
PDB_FILE_NAME = "5ht1a_aripiprazol_model.pdb"

# Insert the name of the PREPI file in between the ""
PREPI_FILE_NAME = "aripiprazol_B3LYP_opt_HF631Gs_charge_nosym.prepi"

# DON'T TOUCH ANYTHING BELOW THIS COMMENT!

pdb_reordered = []
pdb_formatted = []

# Open PDB and create list for ligand
with open('%s' % PDB_FILE_NAME) as pdbfile:
    old_pdb_list = []
    pdb_append_list = []

    for line in pdbfile:
        old_pdb_list.append(line.split())
        if 'HETATM' or 'ATOM' in line:
            if 'MOL' in line:
                pdb_append_list.append(line.split())

# Open MOL2 and create list for ligand
with open('%s' % PREPI_FILE_NAME) as prepifile:
    prepifile_reference_list = []

    # Ignores DUMMY atoms in PREPI
    for line in prepifile:
        if 'DUMM' not in line:
            prepifile_reference_list.append(line.split())

    # Keyword to begin range of PREPI file containing parameter data
    CORRECT_index = [(i, el.index('CORRECT')) for i, el in enumerate(prepifile_reference_list) if 'CORRECT' in el]
    del prepifile_reference_list[:(CORRECT_index[0][0]) + 2]

    # Keyword to end range of PREPI file containing parameter data
    LOOP_index = [(i, el.index('LOOP')) for i, el in enumerate(prepifile_reference_list) if 'LOOP' in el]
    del prepifile_reference_list[LOOP_index[0][0]:]

    # Remove empty elements in nested list
    prepifile_reference_list = [x for x in prepifile_reference_list if x]

# Matching PDB using PREPI as reference
for i in prepifile_reference_list:
    for j in pdb_append_list:
        if i[1] == j[2]:
            pdb_reordered.append(j)

# Appends new atom order from PREPI to old PDB
pdb_index_to_append = [(i, el.index("MOL")) for i, el in enumerate(old_pdb_list) if "MOL" in el]
res_list = [x[0] for x in pdb_index_to_append]

# Dictionary look-up to replace old PDB with reordered PDB based on index from list
replacement_dictionary = dict(zip(res_list, pdb_reordered))
for key, value in replacement_dictionary.items():
    old_pdb_list[key] = value

# Corrects formatting for new PDB file
for i in old_pdb_list:
    print(i)
    pdb_line = []
    if 'TER' in i[0]:
        pdb_line.append(i[0])
        pdb_formatted.append(pdb_line)
    else:
        if 'HETATM' in i[0]:
            pdb_line.append(i[0] + " ")
        if 'ATOM' in i[0]:
            pdb_line.append(i[0] + "   ")
        pdb_line.append(" " * (4 - len(i[1])) + i[1])
        pdb_line.append("  " + i[2])
        pdb_line.append(" " * (4 - len(i[2])) + i[3])
        pdb_line.append(" " + i[4])
        pdb_line.append(" " * (2 - len(i[4])) + i[5])
        pdb_line.append((" " * (8 - len(i[5]))) + i[6])
        pdb_line.append((" " * (8 - len(i[6]))) + i[7])
        pdb_line.append((" " * (8 - len(i[7]))) + i[8])
        pdb_line.append((" " * (9 - len(i[8]))) + i[9])
        pdb_line.append((" " * (5 - len(i[9]))) + i[10])
        pdb_line.append((" " * (16 - len(i[10]))) + i[11])
        pdb_formatted.append(pdb_line)

# Writes a new PDB file with correct ordering
with open('%s_reordered.pdb' % PDB_FILE_NAME, 'w') as f:
    for i in pdb_formatted:
        f.write(str("".join(i)) + "\n")
