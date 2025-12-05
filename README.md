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

## 🚀 Deploy to Vercel

[![Deploy with Vercel](https://vercel.com/button)](https://vercel.com/new/clone?repository-url=https://github.com/khushibaweja/motif-domain)

### Quick Deploy Steps:

1. Click the button above OR go to [vercel.com](https://vercel.com)
2. Import this GitHub repository
3. Vercel will auto-detect Python and deploy
4. Your app will be live at `https://motif-domain.vercel.app`

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
├── vercel.json         # Vercel configuration
└── README.md
```

## 🛠️ Technologies

- **Backend**: Flask, BioPython
- **Visualization**: Matplotlib, Plotly
- **Frontend**: Bootstrap 5
- **API**: NCBI Entrez E-utilities
- **Deployment**: Vercel

## 👤 Author

**Khushi Baweja** - [@khushibaweja](https://github.com/khushibaweja)

---

Made with ❤️ for bioinformatics
