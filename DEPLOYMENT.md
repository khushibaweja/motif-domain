# Deployment Guide - Protein Domain Visualization Tool

## 🚀 Deployment Options

### Option 1: Deploy to Render (Recommended - Free Tier Available)

1. **Prepare your repository**

   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   ```

2. **Push to GitHub**
   - Create a new repository on GitHub
   - Push your code:

   ```bash
   git remote add origin https://github.com/yourusername/protein-viz-tool.git
   git branch -M main
   git push -u origin main
   ```

3. **Deploy on Render**
   - Go to [render.com](https://render.com)
   - Sign up/Login with GitHub
   - Click "New +" → "Web Service"
   - Connect your GitHub repository
   - Configure:
     - **Name**: protein-viz-tool
     - **Environment**: Python 3
     - **Build Command**: `pip install -r requirements-web.txt`
     - **Start Command**: `gunicorn app:app`
   - Add Environment Variable:
     - **Key**: `NCBI_EMAIL`
     - **Value**: your.email@example.com
   - Click "Create Web Service"

Your app will be live at: `https://protein-viz-tool.onrender.com`

---

### Option 2: Deploy to Railway

1. **Push code to GitHub** (same as above)

2. **Deploy on Railway**
   - Go to [railway.app](https://railway.app)
   - Sign up/Login with GitHub
   - Click "New Project" → "Deploy from GitHub repo"
   - Select your repository
   - Railway auto-detects Python and deploys
   - Add Environment Variable:
     - Go to Variables tab
     - Add `NCBI_EMAIL` = your.email@example.com

---

### Option 3: Deploy to Heroku

1. **Install Heroku CLI**

   ```bash
   # Windows (using Chocolatey)
   choco install heroku-cli

   # Or download from: https://devcenter.heroku.com/articles/heroku-cli
   ```

2. **Create Procfile**
   Already created in your project: `Procfile`

3. **Deploy**
   ```bash
   heroku login
   heroku create protein-viz-tool
   heroku config:set NCBI_EMAIL=your.email@example.com
   git push heroku main
   heroku open
   ```

---

### Option 4: Deploy to PythonAnywhere

1. **Sign up** at [pythonanywhere.com](https://www.pythonanywhere.com)

2. **Upload your code**
   - Go to "Files" tab
   - Upload all files or use Git

3. **Install dependencies**
   - Go to "Consoles" → Start a Bash console

   ```bash
   pip install --user -r requirements-web.txt
   ```

4. **Configure Web App**
   - Go to "Web" tab
   - Click "Add a new web app"
   - Choose "Manual configuration" → Python 3.10
   - Set source code directory
   - Edit WSGI file to point to your app
   - Set environment variable in WSGI file:
   ```python
   import os
   os.environ['NCBI_EMAIL'] = 'your.email@example.com'
   ```

   - Reload web app

---

### Option 5: Deploy to Azure App Service

1. **Install Azure CLI**

   ```powershell
   winget install Microsoft.AzureCLI
   ```

2. **Login and create resources**

   ```bash
   az login
   az group create --name protein-viz-rg --location eastus
   az appservice plan create --name protein-viz-plan --resource-group protein-viz-rg --sku F1 --is-linux
   az webapp create --resource-group protein-viz-rg --plan protein-viz-plan --name protein-viz-tool --runtime "PYTHON:3.10"
   ```

3. **Configure and deploy**
   ```bash
   az webapp config appsettings set --resource-group protein-viz-rg --name protein-viz-tool --settings NCBI_EMAIL=your.email@example.com
   az webapp up --name protein-viz-tool --resource-group protein-viz-rg
   ```

---

### Option 6: Deploy to Google Cloud Run

1. **Install Google Cloud SDK**
   Download from: https://cloud.google.com/sdk/docs/install

2. **Create Dockerfile** (already created)

3. **Deploy**
   ```bash
   gcloud auth login
   gcloud config set project YOUR_PROJECT_ID
   gcloud run deploy protein-viz-tool --source . --platform managed --region us-central1 --allow-unauthenticated --set-env-vars NCBI_EMAIL=your.email@example.com
   ```

---

### Option 7: Deploy to Vercel (with Serverless Functions)

1. **Install Vercel CLI**

   ```bash
   npm install -g vercel
   ```

2. **Create vercel.json** (already created)

3. **Deploy**
   ```bash
   vercel login
   vercel --prod
   ```

---

## 🔧 Local Testing

Before deploying, test locally:

```bash
# Install dependencies
pip install -r requirements-web.txt

# Set environment variable
$env:NCBI_EMAIL="your.email@example.com"

# Run the app
python app.py
```

Visit: http://localhost:5000

---

## 📋 Environment Variables

Required for all deployments:

- `NCBI_EMAIL`: Your email for NCBI API compliance

Optional:

- `PORT`: Port number (auto-configured by most platforms)
- `FLASK_ENV`: Set to `production` for production deployments

---

## 🔒 Security Considerations

1. **Never commit sensitive data** - Use environment variables
2. **Enable HTTPS** - Most platforms provide this by default
3. **Set rate limiting** - Protect against API abuse
4. **Monitor usage** - Set up logging and monitoring

---

## 📊 Monitoring

Most platforms provide built-in monitoring. Check:

- Response times
- Error rates
- Memory usage
- API call limits (NCBI has rate limits)

---

## 🆘 Troubleshooting

**Issue**: Module not found

- **Solution**: Ensure `requirements-web.txt` is installed

**Issue**: NCBI API fails

- **Solution**: Check if `NCBI_EMAIL` environment variable is set

**Issue**: Timeout errors

- **Solution**: Increase timeout settings in platform config

**Issue**: Memory errors

- **Solution**: Upgrade to paid tier with more memory

---

## 📞 Support

For issues specific to:

- **Code**: Check GitHub issues or create a new one
- **Deployment**: Refer to platform-specific documentation
- **NCBI API**: Check [NCBI E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25500/)

---

## 🎯 Recommended Platform Comparison

| Platform       | Free Tier  | Ease of Use | Best For                |
| -------------- | ---------- | ----------- | ----------------------- |
| Render         | ✅ Yes     | ⭐⭐⭐⭐⭐  | Beginners, Quick Deploy |
| Railway        | ✅ Yes     | ⭐⭐⭐⭐⭐  | Modern UI, Fast         |
| Heroku         | ⚠️ Limited | ⭐⭐⭐⭐    | Established Platform    |
| PythonAnywhere | ✅ Yes     | ⭐⭐⭐      | Python-specific         |
| Azure          | ⚠️ Credits | ⭐⭐⭐      | Enterprise              |
| Google Cloud   | ⚠️ Credits | ⭐⭐⭐      | Scalability             |
| Vercel         | ✅ Yes     | ⭐⭐⭐⭐    | Serverless              |

**Recommendation**: Start with **Render** or **Railway** for easiest deployment with generous free tiers.
