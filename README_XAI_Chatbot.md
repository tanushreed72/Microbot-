# MicrobiomeBot - XAI Chatbot for Microbiome Analysis

An Explainable AI chatbot that provides transparent insights into microbiome interaction analysis using the InfIntE framework.

## 🧬 Overview

MicrobiomeBot is an interactive web-based chatbot that explains microbiome ecological interactions inferred through abductive logic programming. Built on top of the InfIntE R package, it provides transparent, interpretable explanations of complex microbial relationships.

## ✨ Key Features

### Explainable AI Capabilities
- **Logical Reasoning Explanations**: Understand how interactions are inferred through abductive logic
- **Compression Value Interpretations**: Learn what statistical support means for each interaction
- **Step-by-Step Process Breakdown**: See the complete analysis pipeline from data to insights
- **Natural Language Queries**: Ask questions in plain English about your results

### Interaction Analysis
- **5 Interaction Types**: Competition, Mutualism, Predation, Amensalism, Commensalism
- **Network Visualization**: Interactive plots showing microbial interaction networks
- **Statistical Support**: Quantitative measures backing each inference
- **Bidirectional Relationships**: Understand complex ecological dynamics

### User Interface
- **Modern Web Interface**: Clean, responsive design for optimal user experience
- **Real-time Chat**: Instant responses to your questions
- **Data Upload**: Easy CSV/TSV file upload for analysis
- **Interactive Visualizations**: Hover-enabled network plots with detailed information

## 🚀 Quick Start

### Prerequisites

1. **R Environment**: Ensure R is installed with the InfIntE package
```r
library(devtools)
install_github("didacb/InfIntE")
library(InfIntE)
load_PyGol()  # Compile PyGol for abduction
```

2. **Python Dependencies**: Install required packages
```bash
pip install -r requirements.txt
```

### Installation

1. Clone or download the repository
2. Navigate to the project directory
3. Install Python dependencies:
```bash
pip install -r requirements.txt
```

4. Run the chatbot:
```bash
python microbiome_xai_chatbot.py
```

5. Open your browser to `http://localhost:5000`

## 💬 How to Use

### Starting a Conversation
- Type questions in natural language
- Use the example queries in the sidebar
- Ask about interaction types, analysis processes, or specific relationships

### Uploading Data
1. Click "Upload Data" in the sidebar
2. Select a CSV/TSV file with OTU/ASV abundance data
3. Wait for analysis to complete
4. Ask questions about your results

### Example Queries
- "What are the interaction types?"
- "How does the analysis work?"
- "Explain the interaction between Bacillus and Pseudomonas"
- "Show me the network"
- "What does a compression value of 150 mean?"

## 🔬 Scientific Background

### InfIntE Framework
The chatbot is built on the InfIntE (Inference of Interactions using Explainable machine learning) framework, which:

1. **Converts abundance data to logical clauses**
2. **Applies abductive reasoning via PyGol**
3. **Classifies interactions based on effect patterns**
4. **Uses StARS for robust model selection**

### Interaction Types Explained

| Type | Description | Example |
|------|-------------|---------|
| **Competition** | Both species negatively affect each other | Resource competition |
| **Mutualism** | Both species positively affect each other | Symbiotic relationships |
| **Predation** | One benefits, other is harmed | Predator-prey dynamics |
| **Amensalism** | One is harmed, other unaffected | Antibiotic production |
| **Commensalism** | One benefits, other unaffected | Metabolic byproduct utilization |

### Explainability Features

#### Logical Reasoning
Every interaction inference is backed by:
- **Logical clauses** derived from abundance patterns
- **Abductive explanations** finding most parsimonious relationships
- **Hypothesis rules** defining interaction patterns

#### Statistical Support
- **Compression values** measure explanatory power
- **Higher values** indicate stronger statistical support
- **Model selection** ensures robust interaction identification

## 🛠️ Technical Architecture

### Backend Components
- **Flask Web Server**: Handles HTTP requests and responses
- **R Integration**: Via rpy2 for seamless InfIntE package access
- **Data Processing**: Pandas for data manipulation
- **Network Analysis**: NetworkX for graph operations
- **Visualization**: Plotly for interactive plots

### Frontend Components
- **Responsive Design**: Modern CSS with mobile support
- **Real-time Chat**: JavaScript-based messaging system
- **File Upload**: Client-side CSV/TSV parsing
- **Interactive Plots**: Plotly.js integration

## 📊 Data Format

### Input Requirements
- **CSV or TSV format**
- **Rows**: OTUs/ASVs (taxa)
- **Columns**: Samples
- **Values**: Abundance counts
- **Optional**: Sequencing depth information

### Example Data Structure
```
OTU_ID,Sample1,Sample2,Sample3,Sample4
OTU_001,150,200,180,220
OTU_002,50,80,60,90
OTU_003,300,250,280,240
```

## 🔧 Configuration

### Environment Variables
- `FLASK_ENV`: Set to 'development' for debugging
- `FLASK_PORT`: Custom port (default: 5000)
- `R_HOME`: Path to R installation (if needed)

### Customization Options
- **Interaction thresholds**: Modify in R analysis parameters
- **UI styling**: Edit CSS in `templates/index.html`
- **Response templates**: Customize in chatbot class methods

## 🐛 Troubleshooting

### Common Issues

1. **InfIntE package not found**
   - Ensure R is properly installed
   - Install InfIntE package using devtools
   - Run `load_PyGol()` to compile dependencies

2. **rpy2 connection errors**
   - Check R installation path
   - Verify R packages are accessible
   - Set R_HOME environment variable if needed

3. **File upload failures**
   - Check file format (CSV/TSV)
   - Ensure proper column structure
   - Verify data contains numeric values

4. **Network visualization issues**
   - Ensure analysis has been run first
   - Check browser console for JavaScript errors
   - Verify Plotly.js is loading correctly

## 📚 References

> Barroso-Bergada D, Tamaddoni-Nezhad A, Varghese D, Vacher C, Galic N, Laval V, Suffert F, Bohan DA (2023). "Unravelling the web of dark interactions: Explainable inference of the diversity of microbial interactions." In *Advances in Ecological Research: Roadmaps: Part A*, volume 68, 155-183. Academic Press.

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## 📄 License

This project builds upon the InfIntE package. Please refer to the original package license for usage terms.

## 🆘 Support

For issues or questions:
1. Check the troubleshooting section
2. Review the InfIntE package documentation
3. Open an issue on the repository
4. Contact the development team

---

**Made with 🧬 for the microbiome research community**
