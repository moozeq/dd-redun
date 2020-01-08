#!/usr/bin/env bash
#   $1  -   smi database
#   $2  -   ligand index
#   $3  -   receptor index

echo "[*] smiles database = $1, ligand index = $2, receptor index = $3"

LIGAND_INDEX=$((${2} + 1))
LIGAND_NAME_TEMP=$(sed "${LIGAND_INDEX}q;d" ${1} | awk '{print $2}')
LIGAND_SMI=$(sed "${LIGAND_INDEX}q;d" ${1} | awk '{print $1}')
LIGAND_NAME=${LIGAND_NAME_TEMP:0:${#LIGAND_NAME_TEMP}-7}
LIGAND_FILE_SMI=coreset/${LIGAND_NAME}/${LIGAND_NAME}_ligand.smi
LIGAND_FILE_PDBQT=coreset/${LIGAND_NAME}/${LIGAND_NAME}_ligand.pdbqt

if [[ ! -f "${LIGAND_FILE_SMI}" ]]; then
    echo ${LIGAND_SMI} > ${LIGAND_FILE_SMI}
fi

if [[ ! -f "${LIGAND_FILE_PDBQT}" ]]; then
    obabel -ismi ${LIGAND_FILE_SMI} -opdbqt --gen3d -O ${LIGAND_FILE_PDBQT}
    echo "[+] Converted ligand: smi -> pdbqt"
fi

RECEPTOR_INDEX=$((${3} + 1))
RECEPTOR_NAME_TEMP=$(sed "${RECEPTOR_INDEX}q;d" ${1} | awk '{print $2}')
RECEPTOR_NAME=${RECEPTOR_NAME_TEMP:0:${#RECEPTOR_NAME_TEMP}-7}
RECEPTOR_FILE_PDB=coreset/${RECEPTOR_NAME}/${RECEPTOR_NAME}_protein.pdb
RECEPTOR_FILE_PDBQT=coreset/${RECEPTOR_NAME}/${RECEPTOR_NAME}_protein.pdbqt
RECEPTOR_LIGAND_FILE_MOL2=coreset/${RECEPTOR_NAME}/${RECEPTOR_NAME}_ligand.mol2

if [[ ! -f "${RECEPTOR_FILE_PDBQT}" ]]; then
    obabel -xc -xr -ipdb ${RECEPTOR_FILE_PDB} -opdbqt -O ${RECEPTOR_FILE_PDBQT}
    echo "[+] Converted protein: pdb -> pdbqt"
fi

DOCKING_RESULTS_FILE=docking_ligand${2}_protein${3}.pdbqt
DOCKING_RESULTS_FILE_RESCORED=docking_ligand${2}_protein${3}_rescored.mol2

if [[ ! -f "${DOCKING_RESULTS_FILE}" ]]; then
    COORDS_NUM=$(sed -n -e '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/ p' ${RECEPTOR_LIGAND_FILE_MOL2} | wc -l)
    ((COORDS_NUM = COORDS_NUM - 2)) # remove first and last empty lines
    COORD_X=$(sed -n -e '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/ p' ${RECEPTOR_LIGAND_FILE_MOL2} | awk '{print $3}' | awk '{s+=$1} END {print s}')
    COORD_X=$(echo "$COORD_X / $COORDS_NUM" | bc -l)
    COORD_Y=$(sed -n -e '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/ p' ${RECEPTOR_LIGAND_FILE_MOL2} | awk '{print $4}' | awk '{s+=$1} END {print s}')
    COORD_Y=$(echo "$COORD_Y / $COORDS_NUM" | bc -l)
    COORD_Z=$(sed -n -e '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/ p' ${RECEPTOR_LIGAND_FILE_MOL2} | awk '{print $5}' | awk '{s+=$1} END {print s}')
    COORD_Z=$(echo "$COORD_Z / $COORDS_NUM" | bc -l)

    BOX_X=20
    BOX_Y=20
    BOX_Z=20

    echo "[+] Protein active center coords = (${COORD_X}, ${COORD_Y}, ${COORD_Z})"
    echo "[+] Protein acitve center box size = (${BOX_X}, ${BOX_Y}, ${BOX_Z})"

    vina --receptor ${RECEPTOR_FILE_PDBQT} --ligand ${LIGAND_FILE_PDBQT} --size_x ${BOX_X} --size_y ${BOX_Y} --size_z ${BOX_Z} --center_x ${COORD_X} --center_y ${COORD_Y} --center_z ${COORD_Z} --out ${DOCKING_RESULTS_FILE}
    echo "[*] Docking completed, results at: ${DOCKING_RESULTS_FILE}"
fi

if [[ ! -f "${DOCKING_RESULTS_FILE_RESCORED}" ]]; then
    oddt_cli ${DOCKING_RESULTS_FILE} --receptor ${RECEPTOR_FILE_PDBQT} --score rfscore_v1 --score rfscore_v2 --score rfscore_v3 --score nnscore -O ${DOCKING_RESULTS_FILE_RESCORED}
    echo "[*] Rescoring completed, results at: ${DOCKING_RESULTS_FILE_RESCORED}"
fi
exit 0