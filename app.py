import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import os
import tempfile
import base64

def smiles_to_pdb(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb:
        Chem.MolToPDBFile(mol, temp_pdb.name)
        return temp_pdb.name, mol

def lipinski_rule_of_five(mol):
    violations = 0
    parameters = {
        "Number of Hydrogen Bond Donors (HBD)": Descriptors.NumHDonors(mol),
        "Number of Hydrogen Bond Acceptors (HBA)": Descriptors.NumHAcceptors(mol),
        "Molecular Weight (MW)": Descriptors.MolWt(mol),
        "LogP (lipophilicity)": Descriptors.MolLogP(mol)
    }
    st.write("### Lipinski's Rule of Five Parameters:")
    for parameter, value in parameters.items():
        st.text(f"{parameter}: {value}")
    if parameters["Number of Hydrogen Bond Donors (HBD)"] > 5:
        violations += 1
    if parameters["Number of Hydrogen Bond Acceptors (HBA)"] > 10:
        violations += 1
    if parameters["Molecular Weight (MW)"] > 500:
        violations += 1
    if parameters["LogP (lipophilicity)"] > 5:
        violations += 1
    if violations == 0:
        st.success("Lipinski's Rule of Five: No violations")
    else:
        st.error(f"Lipinski's Rule of Five: {violations} violation(s)")
    st.write("###")

def ghose_rule(mol):
    violations = 0
    parameters = {
        "Molecular Weight (MW)": Descriptors.MolWt(mol),
        "LogP (lipophilicity)": Descriptors.MolLogP(mol),
        "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Number of Aromatic Rings": Descriptors.NumAromaticRings(mol)
    }
    st.write("### Ghose's Rule Parameters:")
    for parameter, value in parameters.items():
        st.text(f"{parameter}: {value}")
    if parameters["Molecular Weight (MW)"] < 480 or parameters["Molecular Weight (MW)"] > 500:
        violations += 1
    if parameters["LogP (lipophilicity)"] < 0.4 or parameters["LogP (lipophilicity)"] > 5.6:
        violations += 1
    if parameters["Number of Rotatable Bonds"] > 10:
        violations += 1
    if parameters["Number of Aromatic Rings"] > 2:
        violations += 1
    if violations == 0:
        st.success("Ghose's Rule: No violations")
    else:
        st.error(f"Ghose's Rule: {violations} violation(s)")
    st.write("###")

def veber_rule(mol):
    violations = 0
    parameters = {
        "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Molecular Weight (MW)": Descriptors.MolWt(mol)
    }
    st.write("### Veber's Rule Parameters:")
    for parameter, value in parameters.items():
        st.text(f"{parameter}: {value}")
    if parameters["Number of Rotatable Bonds"] <= 10 and parameters["Molecular Weight (MW)"] <= 500:
        st.success("Veber's Rule: Passed")
    else:
        st.error("Veber's Rule: Failed")
    st.write("###")

def main():
    st.title("SMILES to PDB Converter and 3D Visualizer")
    st.write("""
    This web-based application allows you to convert SMILES notation to PDB format and visualize the 3D structure. It includes Lipinski's Rule of Five, Ghose's Rule, and Veber's Rule analysis. Developed by Dr. Karthick Vasudevan.
    """)

    smiles_input = st.text_input("Enter SMILES notation:")

    if st.button("Convert to PDB"):
        pdb_file, mol = smiles_to_pdb(smiles_input)
        if pdb_file is not None:
            st.success("Conversion successful! PDB file generated.")
            lipinski_rule_of_five(mol)
            ghose_rule(mol)
            veber_rule(mol)
            st.write("### 3D Visualization:")
            st.write("Click the link below to view the 3D visualization:")
            st.markdown(f"[Open 3D Visualization](https://molview.org/?inputFormat=pdb&structureUrl=data%3Achemical%2Fx-pdb%3Bbase64%2C{base64.b64encode(open(pdb_file, 'rb').read()).decode()})", unsafe_allow_html=True)
            st.write("### Download PDB file:")
            with open(pdb_file, "rb") as f:
                pdb_bytes = f.read()
            st.download_button(
                label="Download PDB file",
                data=pdb_bytes,
                file_name="molecule.pdb",
                mime="chemical/x-pdb"
            )
            st.write("###")
        else:
            st.error("Invalid SMILES notation.")

if __name__ == "__main__":
    main()
