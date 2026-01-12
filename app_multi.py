import streamlit as st
import streamlit.components.v1 as components
import os
import subprocess
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

# ---------------------------------------------------------
# Ρυθμίσεις Σελίδας & Λογότυπο
# ---------------------------------------------------------
st.set_page_config(page_title="Ionian University Docking Tool", layout="wide")

# --- ΠΡΟΣΘΗΚΗ ΛΟΓΟΤΥΠΟΥ ΣΤΗΝ SIDEBAR ---
LOGO_PATH = "image_4.jpg"
if os.path.exists(LOGO_PATH):
    st.sidebar.image(LOGO_PATH, use_container_width=True)
else:
    # Αν δεν βρεθεί η εικόνα, δεν κρασάρει, απλά συνεχίζει
    pass
# ---------------------------------------

st.title("Εργαλείο Μοριακής Πρόσδεσης (Docking)")
st.markdown("Υπολογισμός ενέργειας σύνδεσης (Binding Affinity) με χρήση AutoDock Vina.")

# Έλεγχος αν υπάρχει το Vina
VINA_PATH = "vina.exe"
if not os.path.exists(VINA_PATH):
    st.error(" ΣΦΑΛΜΑ: Λείπει το αρχείο `vina.exe` από τον φάκελο της εφαρμογής.")
    st.stop()

# ---------------------------------------------------------
# Sidebar - Ρυθμίσεις
# ---------------------------------------------------------
st.sidebar.header("1. Εισαγωγή Υποδοχέα (Receptor)")
uploaded_receptor = st.sidebar.file_uploader("Αρχείο Πρωτεΐνης (.pdbqt)", type=["pdbqt"])

st.sidebar.header("2. Ρυθμίσεις Κουτιού (Grid Box)")
st.sidebar.info("Εισάγετε εδώ τις συντεταγμένες που βρήκατε στο Chimera.")

# Session State για να κρατάμε τα νούμερα (default values)
if 'center' not in st.session_state:
    st.session_state.center = (0.0, 0.0, 0.0)
if 'size' not in st.session_state:
    st.session_state.size = (0.0, 0.0, 0.0)

# Inputs για τις συντεταγμένες (Χειροκίνητα)
st.sidebar.subheader("Κέντρο (Center)")
c_x = st.sidebar.number_input("Center X", value=st.session_state.center[0], format="%.3f")
c_y = st.sidebar.number_input("Center Y", value=st.session_state.center[1], format="%.3f")
c_z = st.sidebar.number_input("Center Z", value=st.session_state.center[2], format="%.3f")

st.sidebar.subheader("Μέγεθος (Size)")
s_x = st.sidebar.number_input("Size X", value=st.session_state.size[0], format="%.3f")
s_y = st.sidebar.number_input("Size Y", value=st.session_state.size[1], format="%.3f")
s_z = st.sidebar.number_input("Size Z", value=st.session_state.size[2], format="%.3f")

st.sidebar.subheader("Παράμετροι Vina")
exhaustiveness = st.sidebar.slider("Exhaustiveness", 1, 32, 8,
                                   help="Υψηλότερες τιμές αυξάνουν την ακρίβεια αλλά και τον χρόνο.")

# ---------------------------------------------------------
# Κυρίως Εφαρμογή
# ---------------------------------------------------------
st.subheader("3. Εισαγωγή Φαρμάκου (Ligand)")
smiles = st.text_input("SMILES Code:", value="", placeholder="π.χ. COC1=C(C=C2C(=C1)CC(C2=O)CC3CCN(CC3)CC4=CC=CC=C4)OC")
st.caption("Εισάγετε τον κωδικό SMILES του μορίου.")

run_btn = st.button(" Έναρξη Docking")

if run_btn and uploaded_receptor and smiles:

    # Καθαρισμός παλιών αρχείων
    for f in ["receptor.pdbqt", "ligand.pdbqt", "output.pdbqt"]:
        if os.path.exists(f): os.remove(f)

    # Καθαρισμός Receptor (αφαιρούμε γραμμές CONECT για να μην κολλάει το Vina)
    raw_content = uploaded_receptor.getvalue().decode("utf-8")
    clean_lines = [line for line in raw_content.splitlines() if not line.startswith("CONECT")]
    with open("receptor.pdbqt", "w") as f:
        f.write("\n".join(clean_lines))

    # Προετοιμασία Ligand από SMILES
    with st.spinner("Προετοιμασία δομής φαρμάκου..."):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol: raise ValueError("Μη έγκυρο SMILES string.")
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            prep = MoleculePreparation()
            prep.prepare(mol)
            prep.write_pdbqt_file("ligand.pdbqt")
        except Exception as e:
            st.error(f"Σφάλμα στο SMILES: {e}")
            st.stop()

    # Εκτέλεση Vina
    with st.spinner("Το Docking εκτελείται..."):
        cmd = [
            VINA_PATH,
            "--receptor", "receptor.pdbqt",
            "--ligand", "ligand.pdbqt",
            "--center_x", str(c_x), "--center_y", str(c_y), "--center_z", str(c_z),
            "--size_x", str(s_x), "--size_y", str(s_y), "--size_z", str(s_z),
            "--exhaustiveness", str(exhaustiveness),
            "--out", "output.pdbqt"
        ]

        process = subprocess.run(cmd, capture_output=True, text=True)

        if process.returncode == 0:
            st.success("Το Docking ολοκληρώθηκε!")

            # Εξαγωγή του Score
            affinity = "N/A"
            for line in process.stdout.splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "1":
                    try:
                        affinity = str(float(parts[1]))
                        break
                    except:
                        continue

            st.metric(" Binding Affinity", f"{affinity} kcal/mol")

            # Οπτικοποίηση
            st.subheader("Αποτελέσματα")
            if os.path.exists("output.pdbqt"):
                with open("output.pdbqt", "r") as f:
                    docked_data = f.read()

                col1, col2 = st.columns([3, 1])
                with col1:
                    view = py3Dmol.view(width=700, height=500)
                    with open("receptor.pdbqt", "r") as f:
                        view.addModel(f.read(), "pdbqt")
                    view.setStyle({'model': -1}, {"cartoon": {"color": "gray", "opacity": 0.7}})
                    view.addModel(docked_data, "pdbqt")
                    view.setStyle({'model': -1}, {"stick": {"colorscheme": "greenCarbon"}})
                    view.zoomTo()
                    components.html(view._make_html(), height=500)

                with col2:
                    st.write("Λήψη")
                    st.download_button("Docked Ligand", docked_data, "docked_ligand.pdbqt")
                    with open("receptor.pdbqt", "r") as f: rec_txt = f.read()
                    st.download_button("Full Complex", rec_txt + "\n" + docked_data, "complex.pdbqt")
        else:
            st.error(" Το Docking απέτυχε. Ελέγξτε τις συντεταγμένες του κουτιού.")
            st.code(process.stdout)

elif run_btn:
    st.warning(" Ανεβάστε αρχείο PDBQT και συμπληρώστε το SMILES.")
