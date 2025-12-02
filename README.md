# Protein Domain Visualization Tool

## Overview

This tool analyzes protein sequences and visualizes domains/motifs with both static and interactive plots. Available as both a **command-line tool** and a **web application**.

## Features

- 📁 Upload FASTA files (protein or nucleotide)
- 🔍 Fetch sequences by NCBI accession number
- ✍️ Manual sequence entry
- 📊 Generate static and interactive visualizations
- 🧬 Basic protein feature analysis
- 🌐 Modern web interface with Bootstrap UI
- 📱 Responsive design for mobile devices

## Quick Start

### Web Application (Recommended)

1. Install dependencies:

```bash
pip install -r requirements-web.txt
```

2. Set environment variable:

```bash
# Windows PowerShell
$env:NCBI_EMAIL="your.email@example.com"

# Linux/Mac
export NCBI_EMAIL="your.email@example.com"
```

3. Run the web server:

```bash
python app.py
```

4. Open your browser to: http://localhost:5000

### Command-Line Tool

### Command-Line Tool

1. Install required packages:

```bash
pip install -r requirements.txt
```

2. Set your NCBI email (optional but recommended):

```bash
# Windows PowerShell
$env:NCBI_EMAIL="your.email@example.com"

# Linux/Mac
export NCBI_EMAIL="your.email@example.com"
```

3. Run the program:

```bash
python file1.py
```

Follow the interactive prompts to analyze sequences.

## Deployment

See [DEPLOYMENT.md](DEPLOYMENT.md) for detailed deployment instructions to:

- Render (Free tier, recommended)
- Railway
- Heroku
- PythonAnywhere
- Azure App Service
- Google Cloud Run
- Vercel

## Important Notes

⚠️ **Mock Data**: Currently uses simulated domain data for demonstration. Real protein domain analysis requires:

- InterProScan API setup
- Local InterProScan installation
- Or integration with other bioinformatics tools (Pfam, SMART, etc.)

## Output Files

- `*.png` - Static domain visualization plots
- `*_interactive.html` - Interactive Plotly visualizations

## Improvements Made

- ✅ Removed hard-coded email addresses
- ✅ Added environment variable support
- ✅ Improved error handling throughout
- ✅ Added visualization options (skip or choose type)
- ✅ Better input validation
- ✅ More informative error messages
- ✅ Enhanced mock data warnings
- ✅ Created modern web application with Flask
- ✅ Responsive Bootstrap UI
- ✅ Interactive Plotly visualizations in browser
- ✅ Multiple deployment options

## Project Structure

```
protein-viz-tool/
├── app.py                      # Flask web application
├── file1.py                    # Command-line tool
├── templates/
│   └── index.html             # Web UI
├── static/
│   └── outputs/               # Generated visualizations
├── uploads/                    # Uploaded FASTA files
├── requirements.txt            # CLI dependencies
├── requirements-web.txt        # Web app dependencies
├── Procfile                    # Heroku deployment
├── Dockerfile                  # Docker/Cloud Run deployment
├── vercel.json                # Vercel deployment
├── DEPLOYMENT.md              # Deployment guide
└── README.md                  # This file
```

## Future Enhancements

- [ ] Real InterProScan API integration
- [ ] Unit tests
- [ ] Logging framework
- [ ] Configuration file support
- [ ] Batch processing mode
- [ ] Export results to JSON/CSV
- [ ] Multiple sequence alignment visualization

## License

MIT License
