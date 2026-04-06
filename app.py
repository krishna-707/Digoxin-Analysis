import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdMolDescriptors
import base64
from io import BytesIO
import py3Dmol
from streamlit.components.v1 import html

# --- 1. PAGE CONFIGURATION ---
st.set_page_config(page_title="Digoxin Analysis", layout="wide", initial_sidebar_state="collapsed")

# --- 2. CUSTOM CSS (Neon Blue Background & Glassmorphism) ---
css = """
<style>
    /* Dark background with minimal neon blue pattern */
    .stApp {
        background-color: #050810;
        background-image: 
            radial-gradient(rgba(0, 243, 255, 0.15) 1px, transparent 1px);
        background-size: 25px 25px;
        color: #e2e8f0;
    }
    
    /* Top Right Text Positioning */
    .top-right-info {
        position: absolute;
        top: 20px;
        right: 20px;
        text-align: right;
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        font-size: 14px;
        color: #e2e8f0;
        z-index: 1000;
        background: rgba(10, 15, 30, 0.6);
        padding: 10px 15px;
        border-radius: 8px;
        border: 1px solid rgba(0, 243, 255, 0.1);
        backdrop-filter: blur(5px);
    }

    /* Glassmorphism Panel Styling */
    .glass-panel {
        background: rgba(10, 15, 30, 0.5);
        backdrop-filter: blur(12px);
        -webkit-backdrop-filter: blur(12px);
        border: 1px solid rgba(0, 243, 255, 0.25);
        box-shadow: 0 0 15px rgba(0, 243, 255, 0.15), inset 0 0 10px rgba(0, 243, 255, 0.05);
        border-radius: 12px;
        padding: 20px;
        margin-bottom: 20px;
    }

    /* Explanation Text Styling */
    .explanation-title {
        color: #00f3ff;
        font-size: 1.2rem;
        margin-bottom: 10px;
        display: flex;
        align-items: center;
        gap: 8px;
    }
    
    .explanation-text {
        color: #cbd5e1;
        font-size: 0.95rem;
        line-height: 1.6;
    }
    
    /* Clean headings */
    h1, h2, h3 {
        color: #ffffff;
    }
    .main-title {
        text-align: center;
        margin-top: 40px;
        margin-bottom: 0px;
        text-shadow: 0 0 15px rgba(0,243,255,0.4);
    }
    .sub-title {
        text-align: center;
        color: #94a3b8;
        font-size: 16px;
        margin-bottom: 40px;
    }
</style>
"""
st.markdown(css, unsafe_allow_html=True)

# --- 3. TOP RIGHT INFORMATION ---
top_right_html = """
<div class="top-right-info">
    <strong>Name:</strong> Krishnan T<br>
    <strong>Class/Sec:</strong> AIML-A<br>
    <strong>Reg No:</strong> <span style="color: #ff4b4b;">RA2511026050029</span>
</div>
"""
st.markdown(top_right_html, unsafe_allow_html=True)

# --- 4. HEADER ---
st.markdown("<h1 class='main-title'>Digoxin Molecular Dashboard</h1>", unsafe_allow_html=True)
st.markdown("<div class='sub-title'>Molecular Stereochemistry and Chiral Center Mapping</div>", unsafe_allow_html=True)

# --- 5. DETAILED EXPLANATION SECTION ---
st.markdown("<h2 style='font-size: 24px; border-bottom: 1px solid rgba(0,243,255,0.2); padding-bottom: 10px;'>📘 Detailed Explanation</h2>", unsafe_allow_html=True)

col_exp1, col_exp2, col_exp3 = st.columns(3)

with col_exp1:
    st.markdown("""
    <div class="glass-panel">
        <div class="explanation-title">🔬 What is Stereochemistry?</div>
        <div class="explanation-text">
            Stereochemistry is the study of the three-dimensional arrangement of atoms in molecules. Even if two molecules have the exact same chemical formula and bonds, their spatial arrangement can differ entirely.<br><br>
            In pharmacology, this spatial arrangement is critical because biological receptors are extremely sensitive to the 3D shape of a drug molecule.
        </div>
    </div>
    """, unsafe_allow_html=True)

with col_exp2:
    st.markdown("""
    <div class="glass-panel">
        <div class="explanation-title">💊 About Digoxin</div>
        <div class="explanation-text">
            Digoxin is a potent cardiac glycoside extracted from the foxglove plant (<i>Digitalis lanata</i>). It is primarily used to treat heart failure and arrhythmias.<br><br>
            <b>Key Points:</b>
            <ul>
                <li>Increases the force of heart contractions.</li>
                <li>Inhibits the Na+/K+ ATPase pump in heart cells.</li>
                <li>Has a very narrow therapeutic window.</li>
                <li>Its specific stereochemistry is absolutely essential for binding to its target receptor.</li>
            </ul>
        </div>
    </div>
    """, unsafe_allow_html=True)

with col_exp3:
    st.markdown("""
    <div class="glass-panel">
        <div class="explanation-title">📐 What is R/S Configuration?</div>
        <div class="explanation-text">
            R/S configuration describes the absolute spatial arrangement of atoms around a specific chiral center using the Cahn-Ingold-Prelog (CIP) priority rules.<br><br>
            <ul>
                <li><b>R (Rectus):</b> Clockwise arrangement of priorities.</li>
                <li><b>S (Sinister):</b> Counter-clockwise arrangement.</li>
            </ul>
            These configurations dictate how the Digoxin molecule folds and interacts with the body.
        </div>
    </div>
    """, unsafe_allow_html=True)

# --- 6. RDKIT MOLECULE PROCESSING ---
# Canonical SMILES for Digoxin
digoxin_smiles = "C[C@@H]1[C@H]([C@H](C[C@@H](O1)O[C@@H]2[C@H](O[C@H](C[C@@H]2O)O[C@@H]3[C@H](O[C@H](C[C@@H]3O)O[C@H]4CC[C@]5([C@@H](C4)CC[C@@H]6[C@@H]5C[C@H]([C@]7([C@@]6(CC[C@@H]7C8=CC(=O)OC8)O)C)O)C)C)C)O)O"
mol = Chem.MolFromSmiles(digoxin_smiles)

# Calculate properties
formula = rdMolDescriptors.CalcMolFormula(mol)
exact_mass = rdMolDescriptors.CalcExactMolWt(mol)
chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

# Generate 2D Structure Image
img = Draw.MolToImage(mol, size=(600, 400), kekulize=True)
buffered = BytesIO()
img.save(buffered, format="PNG")
img_b64 = base64.b64encode(buffered.getvalue()).decode()

# Generate 3D Structure using RDKit AllChem
mol_3d = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_3d, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol_3d)
mol_block = Chem.MolToMolBlock(mol_3d)

# Py3DMol rendering function
def render_3d_molecule(mol_block):
    view = py3Dmol.view(width="100%", height=400)
    view.addModel(mol_block, 'sdf')
    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.4}})
    view.setBackgroundColor('#0a0f1e') # Matches glassmorphism dark theme
    view.zoomTo()
    return view._make_html()

# --- 7. STRUCTURES (2D & 3D SIDE-BY-SIDE) ---
st.markdown("<h2 style='font-size: 24px; border-bottom: 1px solid rgba(0,243,255,0.2); padding-bottom: 10px;'>🧬 Molecular Structures</h2>", unsafe_allow_html=True)

col_2d, col_3d = st.columns(2)

with col_2d:
    st.markdown(f"""
    <div class="glass-panel" style="text-align: center;">
        <h3 style="margin-top: 0; font-size: 18px; color: #00f3ff;">🔗 2D Structure</h3>
        <img src="data:image/png;base64,{img_b64}" style="width: 100%; border-radius: 8px; background-color: white; padding: 10px;">
    </div>
    """, unsafe_allow_html=True)

with col_3d:
    st.markdown("""
    <div class="glass-panel" style="text-align: center; height: 100%;">
        <h3 style="margin-top: 0; font-size: 18px; color: #00f3ff;">🧊 3D Structure (Interactive)</h3>
    </div>
    """, unsafe_allow_html=True)
    # Render 3D component slightly overlapping to fit in the column cleanly
    html(render_3d_molecule(mol_block), height=420)

# --- 8. MOLECULAR ANALYSIS & CHIRAL CENTERS ---
st.markdown("<h2 style='font-size: 24px; border-bottom: 1px solid rgba(0,243,255,0.2); padding-bottom: 10px; margin-top: 20px;'>🧪 Molecular Analysis</h2>", unsafe_allow_html=True)

col_metrics, col_table = st.columns([1, 2])

with col_metrics:
    st.markdown(f"""
    <div class="glass-panel" style="text-align: center; padding-top: 30px; padding-bottom: 30px;">
        <h1 style="color: #ff4b4b; font-size: 4.5rem; margin: 0;">21</h1>
        <p style="color: #94a3b8; margin: 0; font-size: 16px; font-weight: bold;">Stereocenters (Chiral Centers)</p>
    </div>
    <div class="glass-panel" style="padding: 15px;">
        <span style="color: #00f3ff; font-weight: bold;">Chemical Formula:</span> {formula}
    </div>
    <div class="glass-panel" style="padding: 15px;">
        <span style="color: #00f3ff; font-weight: bold;">Exact Mass:</span> {exact_mass:.4f} u
    </div>
    """, unsafe_allow_html=True)

with col_table:
    # Extract table data
    table_data = []
    for idx, config in chiral_centers:
        atom_symbol = mol.GetAtomWithIdx(idx).GetSymbol()
        table_data.append({
            "Atom Index": idx,
            "Element": atom_symbol,
            "Configuration": config
        })
    
    st.markdown("<div class='glass-panel'>", unsafe_allow_html=True)
    st.markdown("<h3 style='margin-top: 0; font-size: 18px; color: #00f3ff;'>📊 Chiral Centers Table</h3>", unsafe_allow_html=True)
    df = pd.DataFrame(table_data)
    st.dataframe(df, use_container_width=True, hide_index=True)
    st.markdown("</div>", unsafe_allow_html=True)

# Footer
st.markdown("""
<div style="text-align: left; color: #64748b; font-size: 12px; margin-top: 40px; border-top: 1px solid rgba(0, 243, 255, 0.2); padding-top: 15px;">
    Chemistry Project Submission | Powered by RDKit & Py3DMol
</div>
""", unsafe_allow_html=True)
