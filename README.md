# 🧬 Protein Domain & Motif Analyzer

A web application for analyzing protein sequences to identify domains, motifs, and physicochemical properties with comprehensive visualizations.

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![Flask](https://img.shields.io/badge/Flask-2.0+-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

## 🌟 Features

- **NCBI Integration**: Fetch protein sequences directly using NCBI accession numbers
- **Domain Detection**: Identify protein domains and motifs using pattern-based analysis
- **Comprehensive Visualization Dashboard**:
  - Protein Info Box (ID, Length, MW, pI, Cys count)
  - Amino Acid Class Distribution (Pie Chart)
  - Physicochemical Property Distribution (Bar Chart)
  - Domain Architecture Diagram
  - Kyte-Doolittle Hydropathy Plot
  - Individual Amino Acid Composition
- **Multiple Input Methods**: NCBI Accession, Manual Sequence, FASTA File Upload

## 🚀 Deployment

This project is deployed on Render: https://motif-domain.onrender.com/

To run locally, follow the steps below.

## 📦 Local Installation

```bash
# Clone the repository
git clone https://github.com/khushibaweja/motif-domain.git
cd motif-domain

# Install dependencies
pip install -r requirements.txt

# Run the application
python app.py
```

Visit `http://127.0.0.1:5000` in your browser.

## 🧪 Test Accession Numbers

| Accession      | Protein              | Source   |
| -------------- | -------------------- | -------- |
| `NP_000509`    | Hemoglobin beta      | NCBI API |
| `NP_001265724` | Insulin              | NCBI API |
| `NP_000537`    | p53 Tumor suppressor | NCBI API |
| `NG_049326`    | Sample protein       | Local DB |
| `NP_000671`    | Adrenergic receptor  | Local DB |

## 📁 Project Structure

```
motif-domain/
├── app.py              # Main Flask application
├── templates/
│   └── index.html      # Web interface
├── static/             # Static files
├── requirements.txt    # Python dependencies
└── README.md
```

## 🛠️ Technologies

- **Backend**: Flask, BioPython
- **Visualization**: Matplotlib, Plotly
- **Frontend**: Bootstrap 5
- **API**: NCBI Entrez E-utilities, EBI InterProScan, RCSB PDB
- **Deployment**: Render

## 👤 Author

**Khushi Baweja** - [@khushibaweja](https://github.com/khushibaweja)

---

Made with ❤️ for bioinformatics

© 2025 Khushi Baweja. All rights reserved.
