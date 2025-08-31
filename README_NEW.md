# MicrobiomeBot - Explainable AI Chatbot for Microbiome Analysis

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org)
[![Flask](https://img.shields.io/badge/Flask-2.3.3-green.svg)](https://flask.palletsprojects.com/)
[![R](https://img.shields.io/badge/R-4.0+-red.svg)](https://r-project.org)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

MicrobiomeBot is an advanced **Explainable AI (XAI) chatbot** designed for microbiome interaction analysis. It combines the power of the **InfIntE R package** with an intuitive web interface to provide transparent, interpretable insights into microbial ecological interactions through natural language conversations.

### Key Features

- 🧬 **Microbiome Interaction Analysis** - Infers 5 types of ecological interactions using abductive logic programming
- 🤖 **Conversational AI Interface** - Natural language queries with real-time chat functionality
- 📊 **Interactive Visualizations** - Network graphs, depth analysis plots, and statistical charts
- 🔍 **Explainable AI** - Step-by-step logical reasoning for each interaction inference
- 📄 **Multi-format Support** - CSV, TSV, TXT, and PDF file processing
- 🔬 **Statistical Validation** - Null hypothesis testing and significance analysis
- 📈 **Data Quality Control** - Automated preprocessing and quality assessment

## Architecture

### Backend (`microbiome_xai_chatbot.py`)
- **Flask Web Application** with 20+ API endpoints
- **R Integration** via subprocess calls to InfIntE package
- **Data Processing Pipeline** with intelligent CSV preprocessing
- **File Management System** for uploads and analysis results
- **Network Analysis** using NetworkX and Plotly
- **PDF Processing** with PyPDF2 and PyMuPDF

### Frontend (`chatbot.html`)
- **Single-Page Application** with responsive design
- **Real-time Chat Interface** with typing indicators
- **Interactive Network Visualization** using Plotly.js
- **File Upload Management** with drag-and-drop support
- **Dynamic Content Rendering** with markdown-style formatting
- **Mobile-Responsive Design** with adaptive layouts

## Installation

### Prerequisites

1. **Python 3.8+** with pip
2. **R 4.0+** with the InfIntE package
3. **Git** for cloning the repository

### Step 1: Clone Repository
```bash
git clone <repository-url>
cd MicrobiomeBot
```

### Step 2: Install Python Dependencies
```bash
cd Microbot-
pip install -r requirements.txt
```

### Step 3: Install R Dependencies
```r
# In R console
install.packages("devtools")
devtools::install_github("your-repo/InfIntE")  # Replace with actual InfIntE repository
library(InfIntE)
load_PyGol()  # Compile PyGol for abduction
```

### Step 4: Start the Application
```bash
python microbiome_xai_chatbot.py
```

The application will be available at `http://localhost:5000`

## Usage

### 1. Data Upload Options

#### **Microbiome Data (CSV/TSV)**
- **Format**: Rows = OTUs/ASVs, Columns = Samples, Values = Abundance counts
- **Single File**: Automatic analysis with interaction inference
- **Multiple Files**: Batch processing with comparison capabilities
- **Depth Files**: Special handling for sequencing depth data (filename must contain "depth")

#### **Research Papers (PDF)**
- **Single/Multiple PDFs**: Automatic text extraction and indexing
- **Search Functionality**: Query PDF content through chat interface
- **Literature Integration**: Reference scientific papers in explanations

### 2. Interaction Types Analyzed
- **Competition**: Both species negatively affect each other
- **Mutualism**: Both species positively affect each other
- **Predation**: One species benefits while the other is harmed
- **Amensalism**: One species is harmed while the other is unaffected
- **Commensalism**: One species benefits while the other is unaffected

### 3. Example Queries
```
"What are the interaction types?"
"Show me the network"
"Explain the interaction between Species A and Species B"
"Generate depth plot"
"Compare selected CSVs"
"Help"
```

### 4. Interactive Features

#### **Network Visualization**
- **Click Edges**: Get detailed explanations of specific interactions
- **Hover Details**: See species names, interaction types, and confidence values
- **Legend**: Color-coded interaction types with counts

#### **Depth Analysis**
- **QC Plots**: Distribution, CDF, per-sample, and group-based visualizations
- **Quality Metrics**: Median, mean, quartiles, and flagged low-depth samples
- **Interactive Charts**: Zoom, pan, and hover for detailed information

## API Endpoints

### Core Analysis
- `POST /analyze` - Single file microbiome analysis
- `POST /analyze_multiple` - Multi-file batch analysis
- `POST /analyze_pairwise` - Two-file comparison analysis

### Visualization
- `GET /network` - Generate interaction network visualization
- `GET /depth` - Create sequencing depth plots
- `POST /explain_detailed` - Detailed interaction explanation

### File Management
- `POST /upload_file` - Upload and process files
- `POST /upload_depth` - Upload depth-specific files
- `GET /get_uploaded_files` - List uploaded files
- `POST /delete_csvs` - Delete selected CSV files

### Comparison & Analysis
- `POST /compare_csvs` - Compare two CSV files
- `POST /preview_csvs` - Preview selected CSV files
- `GET /network_null_test` - Statistical null hypothesis testing

## Data Formats

### OTU Table Format
```csv
OTU_ID,Sample1,Sample2,Sample3,...
OTU_001,150,200,175,...
OTU_002,75,100,50,...
```

### Depth File Format
```csv
SampleID,Depth,Group
Sample1,15000,Control
Sample2,18000,Treatment
```

### Taxonomy File Format
```csv
OTU_ID,Kingdom,Phylum,Class,Order,Family,Genus,Species
OTU_001,Bacteria,Proteobacteria,Gammaproteobacteria,...
```

## Scientific Background

### InfIntE Framework
The chatbot implements the complete InfIntE pipeline:

1. **Data Preprocessing**: Handle missing values, normalize abundance data
2. **Logic Clause Generation**: Convert abundance patterns to logical predicates
3. **Abductive Reasoning**: Use PyGol for hypothesis generation and testing
4. **Interaction Classification**: Categorize relationships into ecological types
5. **Model Selection**: Apply StARS for robust interaction identification
6. **Statistical Validation**: Null hypothesis testing and significance assessment

### Explainability Features

#### **Logical Reasoning Chain**
- **Abundance Patterns**: High/low abundance classifications
- **Co-occurrence Rules**: Presence/absence logical relationships  
- **Effect Inference**: Causal relationship determination
- **Compression Scoring**: Statistical support quantification

#### **Statistical Validation**
- **Compression Values**: Higher values indicate stronger statistical support
- **Null Testing**: Random pair comparison for significance assessment
- **Model Selection**: StARS ensures robust interaction identification
- **Cross-validation**: Multiple subsampling for stability

## Configuration

### Environment Variables
```bash
FLASK_SECRET_KEY=your-secret-key-here
DEPTH_CSV=/path/to/default/depth.csv
```

### R Configuration
The application automatically detects R installation in common Windows paths:
- `C:\Program Files\R\R-*\bin\R.exe`
- `C:\Program Files (x86)\R\R-*\bin\R.exe`

## Performance Optimization

### Dataset Size Handling
- **Small datasets** (<1000 features): 5-minute timeout
- **Large datasets** (1000-5000 features): 10-minute timeout  
- **Ultra-large datasets** (>5000 features): 20-minute timeout

### Memory Management
- Automatic data sparsity detection
- Intelligent feature filtering for sparse datasets
- Optimized R script generation with adaptive parameters

## Troubleshooting

### Common Issues

1. **R Not Found**
   ```
   Solution: Install R and add to PATH, or set R_HOME environment variable
   Download R from: https://cran.r-project.org/bin/windows/base/
   ```

2. **InfIntE Package Error**
   ```
   Solution: Install InfIntE package in R using devtools
   library(devtools)
   install_github("didacb/InfIntE")
   ```

3. **Memory Errors**
   ```
   Solution: Reduce dataset size or increase system memory
   Consider using pre-cleaned data with the R cleaning script
   ```

4. **Timeout Errors**
   ```
   Solution: Use pre-cleaned data or reduce feature count
   Check R script execution in console for debugging
   ```

5. **Network Visualization Issues**
   ```
   Solution: Ensure analysis completed successfully
   Check browser console for JavaScript errors
   Verify Plotly.js is loading correctly
   ```

## Development

### Project Structure
```
MicrobiomeBot/
├── Microbot-/
│   ├── microbiome_xai_chatbot.py    # Main Flask application
│   ├── templates/
│   │   └── chatbot.html             # Frontend interface
│   ├── requirements.txt             # Python dependencies
│   ├── data/                        # Sample datasets
│   ├── uploads/                     # User uploaded files
│   └── cache/                       # Analysis cache
├── R/                               # R analysis scripts
├── extra/                           # Additional utilities
└── README.md                        # This file
```

### Adding New Features

1. **Backend**: Add new routes in `microbiome_xai_chatbot.py`
2. **Frontend**: Update JavaScript functions in `chatbot.html`
3. **Analysis**: Extend R integration for new statistical methods

## Dependencies

### Python Packages (requirements.txt)
```
flask==2.3.3
flask-cors==4.0.0
gunicorn==22.0.0
pandas==1.5.3
numpy==1.24.3
plotly==5.15.0
networkx==3.1
PyPDF2==3.0.1
PyMuPDF==1.23.3
werkzeug==2.3.7
python-dateutil==2.8.2
Jinja2==3.1.2
MarkupSafe==2.1.3
click==8.1.7
itsdangerous==2.1.2
tzlocal==5.2
```

### R Packages
- `InfIntE` - Microbiome interaction inference
- `devtools` - Package installation
- `base R` - Core statistical functions

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit changes (`git commit -am 'Add new feature'`)
4. Push to branch (`git push origin feature/new-feature`)
5. Create Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use MicrobiomeBot in your research, please cite:

```bibtex
@software{microbiomebot2024,
  title={MicrobiomeBot: Explainable AI Chatbot for Microbiome Analysis},
  author={[Your Name]},
  year={2024},
  url={https://github.com/your-repo/MicrobiomeBot}
}
```

## Scientific References

> Barroso-Bergada D, Tamaddoni-Nezhad A, Varghese D, Vacher C, Galic N, Laval V, Suffert F, Bohan DA (2023). "Unravelling the web of dark interactions: Explainable inference of the diversity of microbial interactions." In *Advances in Ecological Research: Roadmaps: Part A*, volume 68, 155-183. Academic Press. https://doi.org/10.1016/bs.aecr.2023.09.005

## Support

For questions, issues, or contributions:
- 🐛 Issues: [GitHub Issues](https://github.com/your-repo/MicrobiomeBot/issues)
- 📖 Documentation: [Wiki](https://github.com/your-repo/MicrobiomeBot/wiki)
- 📧 Email: [your-email@domain.com]

## Acknowledgments

- InfIntE R package developers
- Flask and Plotly.js communities
- Microbiome research community

---

**MicrobiomeBot** - Making microbiome analysis transparent, interpretable, and accessible through conversational AI.
