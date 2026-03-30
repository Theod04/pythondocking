import streamlit as st
import streamlit.components.v1 as components
import os
import subprocess
import py3Dmol
import platform
import stat
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

# ---------------------------------------------------------
# Ρυθμίσεις Σελίδας & Λογότυπο
# ---------------------------------------------------------
st.set_page_config(page_title="Ionian University Docking Tool", layout="wide")

# Προσπάθεια εμφάνισης λογοτύπου
LOGO_PATH = "image_4.jpg"
if os.path.exists(LOGO_PATH):
    st.sidebar.image(LOGO_PATH, use_container_width=True)

st.title("Εργαλείο Μοριακής Πρόσδεσης (Docking)")
st.markdown("Υπολογισμός ενέργειας σύνδεσης (Binding Affinity) με χρήση AutoDock Vina.")

# ---------------------------------------------------------
# Έλεγχος & Ρύθμιση του Vina (Πιο ευέλικτος έλεγχος)
# ---------------------------------------------------------
VINA_PATH = None
# Λίστα με πιθανά ονόματα που μπορεί να πήρε το αρχείο στο GitHub
possible_names = ["vina", "vina.exe", "vina_1.2.5_linux_x86_64"]

for name in possible_names:
    if os.path.exists(name):
        VINA_PATH = "./" + name
        break

if VINA_PATH:
    # Αν ΔΕΝ είμαστε σε Windows (άρα είμαστε στο Streamlit Cloud), δίνουμε δικαιώματα εκτέλεσης
    if platform.system() != "Windows":
        try:
            st.info(f"Ρύθμιση εκτελέσιμου: {VINA_PATH}")
            os.chmod(VINA_PATH, os.stat(VINA_PATH).st_mode | stat.S_IEXEC)
        except Exception as e:
            st.warning(f"Προσοχή στα δικαιώματα: {e}")
else:
    st.error("ΣΦΑΛΜΑ: Δεν βρέθηκε το εκτελέσιμο αρχείο του Vina στο φάκελο.")
    st.info(f"Αρχεία που βλέπει το σύστημα: {os.listdir('.')}")
    st.stop()

# ---------------------------------------------------------
# Sidebar - Ρυθμίσεις
# ---------------------------------------------------------
st.sidebar.header("1. Εισαγωγή Υποδοχέα (Receptor)")
uploaded_receptor = st.sidebar.file_uploader("Αρχείο Πρωτεΐνης (.pdbqt)", type=["pdbqt"])

st.sidebar.header("2. Ρυθμίσεις Κουτιού (Grid Box)")
if 'center' not in st.session_state: st.session_state.center = (0.0, 0.0, 0.0)
if 'size' not in st.session_state: st.session_state.size = (15.0, 15.0, 15.0)

st.sidebar.subheader("Κέντρο (Center)")
c_x = st.sidebar.number_input("Center X", value=st.session_state.center[0], format="%.3f")
c_y = st.sidebar.number_input("Center Y", value=st.session_state.center[1], format="%.3f")
c_z = st.sidebar.number_input("Center Z", value=st.session_state.center[2], format="%.3f")

st.sidebar.subheader("Μέγεθος (Size)")
s_x = st.sidebar.number_input("Size X", value=st.session_state.size[0], format="%.3f")
s_y = st.sidebar.number_input("Size Y", value=st.session_state.size[1], format="%.3f")
s_z = st.sidebar.number_input("Size Z", value=st.session_state.size[2], format="%.3f")

st.sidebar.subheader("Παράμετροι Vina")
exhaustiveness = st.sidebar.slider("Exhaustiveness", 1, 32, 8)

# ---------------------------------------------------------
# Κυρίως Εφαρμογή
# ---------------------------------------------------------
st.subheader("3. Εισαγωγή Φαρμάκου (Ligand)")
smiles = st.text_input("SMILES Code:", placeholder="π.χ. CC(=O)OC1=CC=CC=C1C(=O)O")

run_btn = st.button("Έναρξη Docking")

if run_btn and uploaded_receptor and smiles:
    # Καθαρισμός παλιών αρχείων
    for f in ["receptor.pdbqt", "ligand.pdbqt", "output.pdbqt"]:
        if os.path.exists(f): os.remove(f)

    # Αποθήκευση Receptor
    raw_content = uploaded_receptor.getvalue().decode("utf-8")
    with open("receptor.pdbqt", "w") as f:
        f.write(raw_content)

    # Προετοιμασία Ligand
    with st.spinner("Προετοιμασία δομής φαρμάκου..."):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol: raise ValueError("Μη έγκυρο SMILES.")
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
            
            # Εύρεση Affinity Score
            affinity = "N/A"
            for line in process.stdout.splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "1":
                    affinity = parts[1]
                    break
            
            st.metric("Binding Affinity", f"{affinity} kcal/mol")

            # Οπτικοποίηση
            if os.path.exists("output.pdbqt"):
                with open("output.pdbqt", "r") as f:
                    docked_data = f.read()

                view = py3Dmol.view(width=700, height=500)
                with open("receptor.pdbqt", "r") as f:
                    view.addModel(f.read(), "pdbqt")
                view.setStyle({'model': -1}, {"cartoon": {"color": "gray", "opacity": 0.7}})
                view.addModel(docked_data, "pdbqt")
                view.setStyle({'model': -1}, {"stick": {"colorscheme": "greenCarbon"}})
                view.zoomTo()
                components.html(view._make_html(), height=500)
                
                st.download_button("Λήψη Αποτελέσματος (PDBQT)", docked_data, "docked_ligand.pdbqt")
        else:
            st.error("Το Docking απέτυχε. Ελέγξτε τις ρυθμίσεις του Grid Box.")
            st.code(process.stderr)

elif run_btn:
    st.warning("Παρακαλώ ανεβάστε αρχείο PDBQT και SMILES.")
