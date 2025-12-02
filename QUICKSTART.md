# 🚀 Quick Start Guide

## Web Application is Ready!

Your protein domain visualization tool is now running as a web application!

### ✅ What's Been Created

1. **Flask Web Application** (`app.py`)
   - RESTful API endpoints
   - File upload handling
   - NCBI sequence fetching
   - Real-time analysis

2. **Modern Web Interface** (`templates/index.html`)
   - Beautiful Bootstrap UI
   - Three input methods (Manual, Accession, File Upload)
   - Interactive visualizations with Plotly
   - Static visualizations with Matplotlib
   - Responsive mobile design

3. **Deployment Files**
   - `Procfile` - Heroku deployment
   - `Dockerfile` - Docker/Cloud deployment
   - `vercel.json` - Vercel deployment
   - `requirements-web.txt` - Production dependencies

### 🌐 Access the Web App

The server is currently running at:

- **Local**: http://localhost:5000
- **Network**: http://127.0.0.1:5000

Open your browser and visit the URL above!

### 🎯 How to Use

1. **Choose Input Method**:
   - 📝 **Manual Entry**: Type or paste your sequence
   - 🔍 **NCBI Accession**: Enter accession number (e.g., NP_000508.1)
   - 📁 **Upload File**: Upload a FASTA file

2. **Analyze**:
   - Click the analyze button
   - Wait for results (usually 2-5 seconds)

3. **View Results**:
   - Protein features and statistics
   - Detected domains and motifs
   - Static visualization (PNG)
   - Interactive visualization (hover for details)

### 🛑 Stop the Server

Press `Ctrl+C` in the terminal where the server is running.

### 🚀 Deploy to Production

See [DEPLOYMENT.md](DEPLOYMENT.md) for step-by-step deployment guides to:

#### Free Hosting Options:

1. **Render** ⭐ Recommended
   - Free tier available
   - Automatic HTTPS
   - Auto-deploy from GitHub
   - 5-10 minutes to deploy

2. **Railway**
   - $5 free credit monthly
   - Very fast deployment
   - Modern interface

3. **PythonAnywhere**
   - Free tier for small apps
   - Python-focused

#### Command for Quick Deploy (Render):

```bash
# 1. Push to GitHub
git init
git add .
git commit -m "Initial commit"
git remote add origin https://github.com/yourusername/protein-viz-tool.git
git push -u origin main

# 2. Go to render.com
# 3. Connect GitHub repo
# 4. Deploy automatically
```

### 📝 Environment Variables

Before deploying, set:

```bash
NCBI_EMAIL=your.email@example.com
```

This is required for NCBI API compliance.

### 🎨 Features

✅ **Modern UI** - Bootstrap 5 with custom styling
✅ **Fast Analysis** - Real-time sequence processing
✅ **Interactive Plots** - Hover over domains for details
✅ **Mobile Friendly** - Works on phones and tablets
✅ **Error Handling** - Clear error messages
✅ **File Upload** - Drag & drop FASTA files
✅ **NCBI Integration** - Fetch sequences directly
✅ **Export Ready** - Download visualizations

### 🔧 Customization

Edit these files to customize:

- `templates/index.html` - UI and styling
- `app.py` - Backend logic and API
- `requirements-web.txt` - Add more dependencies

### 📦 Next Steps

1. ✅ Test locally (done!)
2. 🚀 Deploy to cloud (see DEPLOYMENT.md)
3. 🌐 Share your URL
4. 📊 Monitor usage
5. 🎨 Customize design

### 🆘 Need Help?

- **Local issues**: Check terminal for error messages
- **Deployment issues**: See DEPLOYMENT.md
- **Feature requests**: Document in README.md

### 🎉 You're All Set!

Your protein domain visualization tool is now:

- ✅ Running locally
- ✅ Ready for deployment
- ✅ User-friendly web interface
- ✅ Production-ready

**Happy analyzing! 🧬**
