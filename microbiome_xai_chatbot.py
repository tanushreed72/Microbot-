import os
import json
import numpy as np
import pandas as pd
from flask import Flask, render_template, request, jsonify, send_from_directory
import subprocess
import tempfile
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from datetime import datetime
import re
from typing import Dict, List, Tuple, Any
import signal
import sys
import PyPDF2
import fitz  # PyMuPDF for better PDF extraction
from werkzeug.utils import secure_filename

# Enable automatic conversion between pandas and R
# pandas2ri.activate()
# numpy2ri.activate()

app = Flask(__name__)

DEPTH_CSV_DEFAULT = os.environ.get("DEPTH_CSV", "/mnt/data/depth.csv")

class MicrobiomeXAIChatbot:
    def __init__(self):
        """Initialize the XAI chatbot with R integration"""
        self.setup_r_environment()
        self.interaction_types = {
            'Competition': 'Both species negatively affect each other',
            'Mutualism': 'Both species positively affect each other', 
            'Predation': 'One species benefits while the other is harmed',
            'Amensalism': 'One species is harmed while the other is unaffected',
            'Commensalism': 'One species benefits while the other is unaffected'
        }
        self.chat_history = []
        self.current_analysis = None
        self.uploaded_files = {
            'csv_files': {},
            'pdf_files': {}
        }
        self.allowed_extensions = {'csv', 'pdf'}
        self.upload_dir = 'uploads'
        os.makedirs(self.upload_dir, exist_ok=True)
        
    def find_r_executable(self):
        """Find R executable on Windows system"""
        possible_paths = [
            r"C:\Program Files\R\R-*\bin\R.exe",
            r"C:\Program Files\R\R-*\bin\x64\R.exe", 
            r"C:\Program Files (x86)\R\R-*\bin\R.exe",
            r"C:\Program Files (x86)\R\R-*\bin\x64\R.exe",
            "R.exe",  # If in PATH
            "R"       # Unix-style
        ]
        
        import glob
        for pattern in possible_paths:
            if '*' in pattern:
                matches = glob.glob(pattern)
                if matches:
                    # Get the latest version
                    return max(matches)
            else:
                try:
                    result = subprocess.run([pattern, '--version'], 
                                          capture_output=True, text=True, timeout=5)
                    if result.returncode == 0:
                        return pattern
                except:
                    continue
        return None
    
    def setup_r_environment(self):
        """Set up R environment and load InfIntE package"""
        try:
            # Find R executable
            self.r_executable = self.find_r_executable()
            if not self.r_executable:
                print("R not found. Please install R or add it to PATH")
                print("Download R from: https://cran.r-project.org/bin/windows/base/")
                self.infinte_loaded = False
                return
            
            print(f"Found R at: {self.r_executable}")
            
            # Test if R is working
            result = subprocess.run([self.r_executable, '--version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                print("R found but not working properly")
                self.infinte_loaded = False
                return
            
            # Try to load InfIntE package
            try:
                result = subprocess.run([self.r_executable, '--slave', '-e', 'library(InfIntE)'], 
                                      capture_output=True, text=True, timeout=30)
                if result.returncode == 0:
                    self.infinte_loaded = True
                    print("InfIntE package loaded successfully")
                else:
                    print(f"InfIntE package not found. Error: {result.stderr}")
                    print("Install InfIntE in R with: devtools::install_github('your-repo/InfIntE')")
                    self.infinte_loaded = False
            except subprocess.TimeoutExpired:
                print("Timeout loading InfIntE package")
                self.infinte_loaded = False
            except Exception as e:
                print(f"Error loading InfIntE: {e}")
                self.infinte_loaded = False
                
        except Exception as e:
            print(f"Error setting up R environment: {e}")
            self.infinte_loaded = False

    def load_depth_csv(self, csv_path: str = DEPTH_CSV_DEFAULT) -> Tuple[pd.DataFrame, Dict[str, str]]:
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Depth CSV not found at: {csv_path}")

        df = pd.read_csv(csv_path)
        df.columns = [str(c).strip() for c in df.columns]

        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        preferred_depth_names = [
            "depth", "sequencing_depth", "sequence_depth", "reads", "read_count",
            "num_reads", "readcount", "total_reads", "library_size"
        ]

        def best_match_depth():
            lowered = {c.lower(): c for c in df.columns}
            for name in preferred_depth_names:
                if name in lowered and pd.api.types.is_numeric_dtype(df[lowered[name]]):
                    return lowered[name]
            if numeric_cols:
                medians = [(c, float(df[c].median())) for c in numeric_cols]
                medians.sort(key=lambda x: x[1], reverse=True)
                return medians[0][0]
            raise ValueError("No numeric depth-like column found in depth.csv")

        depth_col = best_match_depth()

        preferred_sample_names = ["SampleID", "Sample", "sample", "sample_id", "id", "ID", "Name"]
        sample_col = None
        lowered = {c.lower(): c for c in df.columns}
        for name in preferred_sample_names:
            if name.lower() in lowered and not pd.api.types.is_numeric_dtype(df[lowered[name.lower()]]):
                sample_col = lowered[name.lower()]
                break
        if sample_col is None:
            non_num = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
            if non_num:
                sample_col = non_num[0]
            else:
                sample_col = "_Sample"
                df[sample_col] = [f"S{i+1}" for i in range(len(df))]

        group_candidates = ["group", "condition", "batch", "site", "subject", "cohort", "treatment"]
        group_col = None
        for name in group_candidates:
            if name.lower() in lowered and not pd.api.types.is_numeric_dtype(df[lowered[name.lower()]]):
                group_col = lowered[name.lower()]
                break

        df = df[[c for c in [sample_col, depth_col, group_col] if c is not None]].copy()
        df = df.rename(columns={sample_col: "Sample", depth_col: "Depth"})
        if group_col:
            df = df.rename(columns={group_col: "Group"})

        df["Depth"] = pd.to_numeric(df["Depth"], errors="coerce")
        df = df.dropna(subset=["Depth"])
        df = df[df["Depth"] >= 0]

        schema = {"sample_col": "Sample", "depth_col": "Depth"}
        if "Group" in df.columns:
            schema["group_col"] = "Group"

        self.depth_df = df
        self.depth_schema = schema
        return df, schema

    def _depth_qc_summary(self, df: pd.DataFrame) -> Dict[str, Any]:
        depth = df["Depth"].astype(float)
        q1 = float(np.percentile(depth, 25))
        q3 = float(np.percentile(depth, 75))
        iqr = q3 - q1
        low_cut = q1 - 1.5 * iqr if iqr > 0 else q1 * 0.5
        flagged = df[df["Depth"] < low_cut].copy()

        summary = {
            "n_samples": int(len(df)),
            "min": float(depth.min()),
            "q1": q1,
            "median": float(np.median(depth)),
            "mean": float(depth.mean()),
            "q3": q3,
            "max": float(depth.max()),
            "std": float(depth.std(ddof=1)) if len(depth) > 1 else 0.0,
            "low_depth_threshold": float(low_cut),
            "n_flagged_low_depth": int(len(flagged)),
        }
        if not flagged.empty:
            summary["flagged_samples_preview"] = flagged.sort_values("Depth").head(20).to_dict("records")
        return summary

    def generate_depth_figures(self) -> Dict[str, Any]:
        if self.depth_df is None:
            self.load_depth_csv(DEPTH_CSV_DEFAULT)

        df = self.depth_df.copy()
        qc = self._depth_qc_summary(df)

        hist = px.histogram(
            df, x="Depth", nbins=min(60, max(10, int(np.sqrt(len(df))))),
            opacity=0.9, marginal="box", title="Sequencing Depth — Distribution"
        )
        hist.update_layout(
            bargap=0.02, xaxis_title="Depth (reads)", yaxis_title="Number of samples",
            hovermode="x unified", template="plotly_white", title_x=0.5
        )
        hist.add_vline(x=qc["median"], line_width=2, line_dash="dash",
                       annotation_text=f"Median: {int(qc['median'])}")

        bar_df = df.sort_values("Depth", ascending=False).reset_index(drop=True)
        bar = px.bar(bar_df, x="Sample", y="Depth", title="Sequencing Depth — Per Sample (sorted)")
        bar.update_layout(xaxis_title="Sample", yaxis_title="Depth (reads)",
                          template="plotly_white", title_x=0.5, margin=dict(l=40, r=10, t=60, b=120))

        fig_group_json = None
        if "Group" in df.columns:
            box = px.box(df, x="Group", y="Depth", points="all", title="Sequencing Depth — By Group")
            box.update_layout(xaxis_title="Group", yaxis_title="Depth (reads)",
                              template="plotly_white", title_x=0.5)
            fig_group_json = box.to_json()

        sorted_depth = np.sort(df["Depth"].values.astype(float))
        cdf_y = np.arange(1, len(sorted_depth) + 1) / len(sorted_depth)
        cdf_fig = go.Figure()
        cdf_fig.add_trace(go.Scatter(x=sorted_depth, y=cdf_y, mode="lines", name="CDF"))
        cdf_fig.update_layout(title="Sequencing Depth — Empirical CDF",
                              xaxis_title="Depth (reads)", yaxis_title="Fraction of samples ≤ depth",
                              template="plotly_white", title_x=0.5)

        return {
            "fig_distribution": hist.to_json(),
            "fig_per_sample": bar.to_json(),
            "fig_by_group": fig_group_json,
            "fig_cdf": cdf_fig.to_json(),
            "qc_summary": qc,
            "schema": self.depth_schema or {"sample_col": "Sample", "depth_col": "Depth"}
        }

    
    def analyze_microbiome_data(self, otu_data: pd.DataFrame, depth: List[float] = None) -> Dict[str, Any]:
        """Run InfIntE analysis on microbiome data with optimized performance"""
        if not self.infinte_loaded:
            return {"error": "InfIntE package not available"}
        
        try:
            # Preprocess the data to handle missing values and ensure proper format
            processed_data = self.preprocess_otu_data(otu_data)
            if isinstance(processed_data, dict) and 'error' in processed_data:
                return processed_data
            
            # Create R script for analysis with optimized parameters
            input_path = os.path.join(os.getcwd(), 'input.csv')
            processed_data.to_csv(input_path, index=True)
            
            interactions_output = os.path.join(os.getcwd(), 'interactions.csv')
            interactions_enhanced_output = os.path.join(os.getcwd(), 'interactions_enhanced.csv')
            
            # Fix Windows compatibility: mclapply doesn't support multiple cores on Windows
            if os.name == 'nt':  # Windows
                ncores = 1
            else:
                ncores = min(6, os.cpu_count() or 1)  # Unix/Linux systems
            
            nperms = 10  # Reduced from 20 for faster processing
            
            depth_str = ', '.join(map(str, depth)) if depth else ''
            
            # Convert paths to forward slashes for R compatibility
            input_path_r = input_path.replace('\\', '/')
            interactions_output_r = interactions_output.replace('\\', '/')
            interactions_enhanced_output_r = interactions_enhanced_output.replace('\\', '/')
            
            r_script = f"""
            library(InfIntE)
            otu_data <- read.csv("{input_path_r}", row.names=1, check.names=FALSE)
            
            # Set options for Windows compatibility
            options(mc.cores = {ncores})
            
            # Handle Python3 dependency gracefully
            tryCatch({{
                # Try to run analysis with error handling
                result <- infinte(otu_tb=otu_data, exclusion=TRUE, ncores={ncores}, nperms={nperms}{', depth=c(' + depth_str + ')' if depth is not None else ''})
            }}, error = function(e) {{
                # If Python3 error, try with minimal settings
                cat("Warning: Python3 dependency issue detected, using fallback mode\\n")
                result <<- infinte(otu_tb=otu_data, exclusion=TRUE, ncores=1, nperms=5{', depth=c(' + depth_str + ')' if depth is not None else ''})
            }})
            
            interactions <- result$selected_interactions
            write.csv(interactions, "{interactions_output_r}", row.names=FALSE)
            
            # Also save enhanced interactions with metadata
            interactions_enhanced <- data.frame(
                Species1 = interactions$sp1,
                Species2 = interactions$sp2,
                Interaction_Type = interactions$lnk,
                Compression_Value = interactions$comp,
                Interaction_Strength = cut(interactions$comp, 
                                          breaks = quantile(interactions$comp, c(0, 0.33, 0.66, 1)), 
                                          labels = c("Weak", "Medium", "Strong"),
                                          include.lowest = TRUE),
                Direction = ifelse(grepl("up|app", interactions$lnk), "Positive", "Negative"),
                Timestamp = Sys.time()
            )
            
            write.csv(interactions_enhanced, "{interactions_enhanced_output_r}", row.names=FALSE)
            """
            
            # Write R script
            with open('script.R', 'w') as f:
                f.write(r_script)
            
            # Execute R script with reduced timeout for faster processing
            result = subprocess.run([self.r_executable, '--slave', '-f', 'script.R'], 
                                  capture_output=True, text=True, timeout=180, cwd=os.getcwd())  # Reduced from 300s
            
            if result.returncode != 0:
                return {"error": f"R analysis failed: {result.stderr}"}
            
            # Read results from existing interactions.csv file
            interactions_file = os.path.join(os.getcwd(), 'interactions.csv')
            if os.path.exists(interactions_file):
                interactions = pd.read_csv(interactions_file)
                
                # Merge with existing analysis if available
                if self.current_analysis and 'interactions' in self.current_analysis:
                    # Combine with existing interactions
                    existing_interactions = self.current_analysis['interactions']
                    combined_interactions = pd.concat([existing_interactions, interactions], ignore_index=True)
                    # Remove duplicates based on species pairs
                    combined_interactions = combined_interactions.drop_duplicates(subset=['sp1', 'sp2'], keep='last')
                    interactions = combined_interactions
                
                self.current_analysis = {
                    'interactions': interactions,
                    'otu_data': processed_data,
                    'timestamp': datetime.now().isoformat()
                }
                
                return {
                    'success': True,
                    'interactions': interactions.to_dict('records'),
                    'num_interactions': len(interactions),
                    'interaction_types': interactions['lnk'].value_counts().to_dict() if 'lnk' in interactions.columns else {}
                }
            else:
                return {"error": "No results generated by R analysis"}
                    
        except subprocess.TimeoutExpired:
            return {"error": "Analysis timed out after 3 minutes"}
        except Exception as e:
            return {"error": f"Analysis failed: {str(e)}"}
    
    def preprocess_otu_data(self, otu_data: pd.DataFrame) -> pd.DataFrame:
        """Preprocess data to handle missing values and ensure proper format for any CSV file"""
        try:
            # Make a copy to avoid modifying original data
            processed_data = otu_data.copy()
            
            # Check if data is empty
            if processed_data.empty:
                return {"error": "Empty dataset provided"}
            
            # Check if this is already a cleaned OTU table
            if 'otu_table_cleaned.csv' in str(getattr(self, 'current_file', '')):
                print("Using pre-cleaned OTU data - skipping aggressive filtering")
                skip_cleaning = True
            else:
                skip_cleaning = False
            
            # Handle different CSV formats more flexibly
            # If first column looks like sample names/IDs and has unique values, use it as index
            first_col = processed_data.columns[0]
            if (processed_data[first_col].dtype == 'object' and 
                processed_data[first_col].nunique() == len(processed_data) and
                len(processed_data) > 1):
                processed_data = processed_data.set_index(first_col)
            
            # Remove completely empty rows and columns
            processed_data = processed_data.dropna(how='all', axis=0)  # Remove empty rows
            processed_data = processed_data.dropna(how='all', axis=1)  # Remove empty columns
            
            # Check if we still have data after removing empty rows/columns
            if processed_data.empty:
                return {"error": "Dataset contains only empty values"}
            
            # Handle missing values by filling with 0 (common practice for abundance data)
            processed_data = processed_data.fillna(0)
            
            # Try to convert columns to numeric, but be flexible with mixed data types
            numeric_cols = []
            for col in processed_data.columns:
                try:
                    # Try to convert to numeric
                    numeric_data = pd.to_numeric(processed_data[col], errors='coerce')
                    if not numeric_data.isna().all():  # If at least some values are numeric
                        processed_data[col] = numeric_data.fillna(0)
                        numeric_cols.append(col)
                    else:
                        # For non-numeric columns, try to encode them as categories
                        if processed_data[col].dtype == 'object':
                            # Convert categorical data to numeric codes
                            processed_data[col] = pd.Categorical(processed_data[col]).codes
                            numeric_cols.append(col)
                except:
                    # If conversion fails completely, fill with 0
                    processed_data[col] = 0
                    numeric_cols.append(col)
            
            # Keep only columns that could be converted to numeric
            if numeric_cols:
                processed_data = processed_data[numeric_cols]
            
            # Ensure all values are non-negative for analysis
            processed_data = processed_data.abs()
            
            # Round to reasonable precision
            processed_data = processed_data.round(6)
            
            # Only apply aggressive cleaning if not using pre-cleaned data
            if not skip_cleaning:
                print(f"Original data shape: {processed_data.shape}")
                
                # Calculate sparsity
                total_cells = processed_data.size
                zero_cells = (processed_data == 0).sum().sum()
                sparsity = (zero_cells / total_cells) * 100
                print(f"Data sparsity: {sparsity:.1f}% zeros")
                
                # If dataset is very sparse (>90% zeros) and large (>500 features), apply aggressive filtering
                if sparsity > 90 and processed_data.shape[0] > 500:
                    print("Large sparse dataset detected. Consider using the R cleaning script first:")
                    print("Run: Rscript clean_otu_data.R")
                    print("Then upload the generated 'otu_table_cleaned.csv' file instead.")
                    
                    # Apply moderate filtering for now
                    min_samples = max(2, int(processed_data.shape[1] * 0.02))  # 2% threshold
                    feature_presence = (processed_data > 0).sum(axis=1)
                    abundant_features = feature_presence >= min_samples
                    processed_data = processed_data[abundant_features]
                    print(f"Applied basic filtering: kept {abundant_features.sum()} features")
            
            # Remove rows with all zeros (no data)
            processed_data = processed_data.loc[~(processed_data == 0).all(axis=1)]
            
            # Remove columns with all zeros (no variation)
            processed_data = processed_data.loc[:, ~(processed_data == 0).all(axis=0)]
            
            print(f"Final data shape: {processed_data.shape}")
            
            # More flexible minimum data requirements
            if processed_data.shape[0] < 1:
                return {"error": "No valid data rows found after preprocessing"}
            
            if processed_data.shape[1] < 1:
                return {"error": "No valid data columns found after preprocessing"}
            
            # For interaction analysis, we need at least 2 features and 2 samples
            # But be more flexible about the requirements
            if processed_data.shape[0] < 2 and processed_data.shape[1] >= 2:
                # If we have few rows but multiple columns, transpose the data
                # This handles cases where features are in columns instead of rows
                processed_data = processed_data.T
            
            # Very lenient validation - just ensure we have some data structure
            if processed_data.shape[0] < 1:
                return {"error": "Need at least 1 feature (row) for analysis"}
            
            if processed_data.shape[1] < 1:
                return {"error": "Need at least 1 sample (column) for analysis"}
            
            # Even more lenient variation check
            total_non_zero = (processed_data != 0).sum().sum()
            if total_non_zero < 1:  # At least 1 non-zero value
                return {"error": "Dataset contains only zero values"}
            
            # Ensure row and column names are strings (required by R)
            processed_data.index = processed_data.index.astype(str)
            processed_data.columns = processed_data.columns.astype(str)
            
            # Final check: ensure no infinite or extremely large values
            processed_data = processed_data.replace([np.inf, -np.inf], 0)
            processed_data = processed_data.clip(upper=1e6)  # Cap at reasonable maximum
            
            return processed_data
            
        except Exception as e:
            return {"error": f"Data preprocessing failed: {str(e)}"}
    
    def explain_interaction(self, species1: str, species2: str) -> str:
        """Provide detailed explanation of interaction between two species"""
        if not self.current_analysis:
            return "No analysis data available. Please run an analysis first."
        
        interactions = self.current_analysis['interactions']
        
        # Find interaction between species
        interaction = interactions[
            ((interactions['sp1'] == species1) & (interactions['sp2'] == species2)) |
            ((interactions['sp1'] == species2) & (interactions['sp2'] == species1))
        ]
        
        if interaction.empty:
            return f"No direct interaction found between {species1} and {species2}."
        
        interaction_type = interaction.iloc[0]['lnk']
        compression = interaction.iloc[0].get('comp', 'N/A')
        
        explanation = f"""
**Interaction between {species1} and {species2}:**

**Type:** {interaction_type}
**Definition:** {self.interaction_types.get(interaction_type, 'Unknown interaction type')}
**Compression Value:** {compression}

**Logical Explanation:**
This interaction was inferred through abductive reasoning using PyGol. The system:
1. Converted abundance data into logical clauses
2. Applied abduction to find the most parsimonious explanations
3. Classified the interaction based on effect patterns

**Interpretation:**
- Higher compression values indicate stronger statistical support
- The interaction type reflects the ecological relationship pattern observed in the data
- This inference is based on abundance changes across samples
        """
        
        return explanation.strip()
    
    def explain_analysis_process(self) -> str:
        """Explain the InfIntE analysis process step by step"""
        return """
**InfIntE Analysis Process Explanation:**

**1. Data Preparation**
- OTU/ASV abundance table is processed
- Sequencing depth information is incorporated
- Data is transformed for logical inference

**2. Logic Clause Generation**
- Abundance data is converted into logical predicates
- Presence/absence patterns are encoded
- Abundance changes between samples are captured

**3. Abductive Reasoning (PyGol)**
- Bottom clauses are generated from the data
- Hypothesis rules define interaction patterns
- Abduction finds the most parsimonious explanations
- Compression values measure explanatory power

**4. Interaction Classification**
- Effects are classified into interaction types
- Bidirectional relationships are identified
- Final network of ecological interactions is built

**5. Model Selection**
- StARS (Stability Approach to Regularization Selection) is used
- Multiple subsamples ensure robust interaction selection
- Stability threshold determines final interaction set

**Why This Approach is Explainable:**
- Every inference is backed by logical rules
- Compression values provide quantitative support
- The reasoning process is transparent and traceable
- Biological interpretations are grounded in ecological theory
        """
    
    def generate_network_visualization(self) -> Dict[str, Any]:
        """Generate interactive network visualization of interactions"""
        if not self.current_analysis:
            return {"error": "No analysis data available"}
        
        interactions = self.current_analysis['interactions']
        
        # Create NetworkX graph
        G = nx.from_pandas_edgelist(
            interactions, 
            source='sp1', 
            target='sp2', 
            edge_attr=['lnk', 'comp'],
            create_using=nx.DiGraph()
        )
        
        # Get node positions with better spacing
        pos = nx.spring_layout(G, k=2, iterations=100, seed=42)
        
        # Create interaction type color mapping
        interaction_colors = {
            'Competition': '#ff4444',
            'Mutualism': '#44ff44', 
            'Predation': '#ff8844',
            'Amensalism': '#8844ff',
            'Commensalism': '#44ffff',
            'up': '#44ff44',
            'down': '#ff4444',
            'appearance': '#ffff44',
            'disappearance': '#ff44ff'
        }
        
        # Prepare edge data with colors and better hover info
        edge_x = []
        edge_y = []
        edge_colors = []
        edge_info = []
        
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            
            edge_data = G.edges[edge]
            interaction_type = edge_data.get('lnk', 'unknown')
            compression = edge_data.get('comp', 'N/A')
            
            # Get color for interaction type
            color = interaction_colors.get(interaction_type, '#888888')
            edge_colors.extend([color, color, None])
            
            edge_info.append(f"<b>{edge[0]} → {edge[1]}</b><br>Type: {interaction_type}<br>Compression: {compression}")
        
        # Create edge traces by interaction type for legend
        edge_traces = []
        unique_interactions = interactions['lnk'].unique()
        
        for interaction_type in unique_interactions:
            type_edges = interactions[interactions['lnk'] == interaction_type]
            type_x = []
            type_y = []
            
            for _, row in type_edges.iterrows():
                if row['sp1'] in pos and row['sp2'] in pos:
                    x0, y0 = pos[row['sp1']]
                    x1, y1 = pos[row['sp2']]
                    type_x.extend([x0, x1, None])
                    type_y.extend([y0, y1, None])
            
            edge_trace = go.Scatter(
                x=type_x, y=type_y,
                line=dict(width=3, color=interaction_colors.get(interaction_type, '#888888')),
                hoverinfo='none',
                mode='lines',
                name=f'{interaction_type} ({len(type_edges)})',
                showlegend=True
            )
            edge_traces.append(edge_trace)
        
        # Create node trace with better styling
        node_x = []
        node_y = []
        node_text = []
        node_info = []
        node_sizes = []
        
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            
            # Use actual species name (truncate if too long)
            display_name = str(node)
            if len(display_name) > 15:
                display_name = display_name[:12] + "..."
            node_text.append(display_name)
            
            # Calculate node metrics
            in_degree = G.in_degree(node)
            out_degree = G.out_degree(node)
            total_degree = in_degree + out_degree
            
            # Size nodes based on degree
            node_size = max(20, min(50, 20 + total_degree * 3))
            node_sizes.append(node_size)
            
            node_info.append(f"<b>Species: {node}</b><br>Incoming interactions: {in_degree}<br>Outgoing interactions: {out_degree}<br>Total degree: {total_degree}")
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            hoverinfo='text',
            text=node_text,
            textposition="middle center",
            textfont=dict(size=10, color='white'),
            hovertext=node_info,
            marker=dict(
                size=node_sizes,
                color='#2E86AB',
                line=dict(width=2, color='white'),
                opacity=0.8
            ),
            name='Species',
            showlegend=True
        )
        
        # Create figure with better layout
        fig = go.Figure(data=edge_traces + [node_trace],
                       layout=go.Layout(
                           title=dict(
                               text='Microbiome Interaction Network',
                               x=0.5,
                               font=dict(size=20)
                           ),
                           showlegend=True,
                           hovermode='closest',
                           margin=dict(b=20,l=5,r=5,t=60),
                           annotations=[ dict(
                               text="Hover over nodes and edges for details. Node size indicates interaction degree.",
                               showarrow=False,
                               xref="paper", yref="paper",
                               x=0.005, y=-0.002,
                               xanchor="left", yanchor="bottom",
                               font=dict(color="#666", size=12)
                           )],
                           xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                           yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                           plot_bgcolor='white',
                           paper_bgcolor='white',
                           legend=dict(
                               yanchor="top",
                               y=0.99,
                               xanchor="left",
                               x=0.01
                           )
                       ))
        
        return {
            'plot_json': fig.to_json(),
            'network_stats': {
                'nodes': len(G.nodes()),
                'edges': len(G.edges()),
                'density': nx.density(G),
                'avg_clustering': nx.average_clustering(G.to_undirected()),
                'interaction_types': interactions['lnk'].value_counts().to_dict()
            }
        }
    
    def process_user_query(self, query: str) -> str:
        """Process natural language queries about microbiome interactions"""
        query_lower = query.lower()

        # ---- Depth viewer triggers ----
        depth_terms = [
            "depth plot", "show depth", "sequencing depth", "depth distribution",
            "library size", "depth visualization", "depth stats", "depth qc"
        ]
        if any(term in query_lower for term in depth_terms):
            return "SHOW_DEPTH"
  
        
        # Debug: Check if current_analysis exists
        print(f"DEBUG: Query = '{query_lower}'")
        print(f"DEBUG: current_analysis exists = {self.current_analysis is not None}")
        if self.current_analysis:
            print(f"DEBUG: interactions shape = {self.current_analysis['interactions'].shape}")
        
        # Dataset filtering commands
        if any(word in query_lower for word in ['filter', 'show me', 'find', 'top']) and self.current_analysis:
            try:
                if 'abundance' in query_lower:
                    # Extract abundance threshold
                    import re
                    abundance_match = re.search(r'abundance\s*[><=]+\s*([0-9.]+)', query_lower)
                    if abundance_match:
                        threshold = float(abundance_match.group(1))
                        # Filter interactions by compression value (proxy for abundance)
                        interactions = self.current_analysis['interactions']
                        if 'comp' in interactions.columns:
                            filtered = interactions[interactions['comp'] > threshold/100]  # Scale threshold
                            return f"Found {len(filtered)} interactions with compression > {threshold/100:.3f}. Showing top 10:\n" + \
                                   filtered.head(10)[['sp1', 'sp2', 'lnk', 'comp']].to_string(index=False)
                        else:
                            return "Compression data not available for abundance filtering."
                            
                elif 'top' in query_lower and 'interaction' in query_lower:
                    # Extract number for "top N interactions"
                    import re
                    num_match = re.search(r'top\s+(\d+)', query_lower)
                    if num_match:
                        n = int(num_match.group(1))
                        interactions = self.current_analysis['interactions']
                        if 'comp' in interactions.columns:
                            top_interactions = interactions.nlargest(n, 'comp')
                            return f"Top {n} interactions by compression value:\n" + \
                                   top_interactions[['sp1', 'sp2', 'lnk', 'comp']].to_string(index=False)
                        else:
                            return "Compression data not available for ranking."
                    else:
                        return "Please specify a number, e.g., 'top 10 interactions'"
                
                elif 'top' in query_lower:
                    # Handle cases like "top 10 interactions" without the word "interaction"
                    import re
                    num_match = re.search(r'top\s+(\d+)', query_lower)
                    if num_match:
                        n = int(num_match.group(1))
                        interactions = self.current_analysis['interactions']
                        if 'comp' in interactions.columns:
                            top_interactions = interactions.nlargest(n, 'comp')
                            return f"Top {n} interactions by compression value:\n" + \
                                   top_interactions[['sp1', 'sp2', 'lnk', 'comp']].to_string(index=False)
                        else:
                            return "Compression data not available for ranking."
                    else:
                        return "Please specify a number, e.g., 'top 10 interactions'"
            except Exception as e:
                return f"Error processing filter request: {str(e)}. Please try uploading and analyzing data first."
        
        elif any(word in query_lower for word in ['filter', 'show me', 'find']) and not self.current_analysis:
            return "No analysis data available. Please upload and analyze data first."
        
        # Analysis commands
        if any(word in query_lower for word in ['analyze', 'analysis', 'run analysis']):
            if 'interaction' in query_lower:
                return "To analyze interactions, please use the 'Upload & Analyze' button in the web interface to upload your CSV data file. Once uploaded, I can help you filter and explore the results!"
        
        # Query patterns
        if any(word in query_lower for word in ['explain', 'what is', 'how does']):
            if 'process' in query_lower or 'analysis' in query_lower:
                return self.explain_analysis_process()
            elif 'interaction' in query_lower:
                # Extract species names if mentioned
                species_pattern = r'\b([A-Z][a-z]+)\b'
                species = re.findall(species_pattern, query)
                if len(species) >= 2:
                    return self.explain_interaction(species[0], species[1])
                else:
                    return "Please specify two species names to explain their interaction."
        
        elif 'network' in query_lower or 'visualization' in query_lower or 'show me the network' in query_lower:
            if not self.current_analysis:
                return "No analysis data available. Please upload and analyze data first."
            return "SHOW_NETWORK"
        
        elif 'types' in query_lower and 'interaction' in query_lower:
            types_explanation = "**Interaction Types in InfIntE:**\n\n"
            for interaction_type, description in self.interaction_types.items():
                types_explanation += f"**{interaction_type}:** {description}\n\n"
            return types_explanation
        
        elif 'help' in query_lower:
            return """
**Available Commands:**
- Ask about interaction types: "What are the interaction types?"
- Explain specific interactions: "Explain interaction between SpeciesA and SpeciesB"
- Learn about the process: "How does the analysis work?"
- Request visualizations: "Show me the network"
- Get help: "help"

**Example Queries:**
- "What is mutualism?"
- "Explain the interaction between Bacillus and Pseudomonas"
- "How does InfIntE work?"
- "Show me interaction types"
            """
        
        # General educational queries about microbiome analysis
        elif any(word in query_lower for word in ['logical reasoning', 'reasoning', 'inference', 'inferences']):
            return """**Logical Reasoning Behind Microbiome Inferences:**

InfIntE uses **abductive logic programming** to infer microbial interactions:

1. **Hypothesis Formation**: Creates logical rules about how microbes might interact
2. **Evidence Analysis**: Examines abundance patterns across samples
3. **Logical Deduction**: Uses PyGol to find the best explanations
4. **Statistical Validation**: Applies StARS model selection to ensure robust results

The system reasons: "If species A affects species B, then we should see specific abundance patterns when A is present vs absent."
            """
        
        elif any(word in query_lower for word in ['statistical significance', 'significance', 'statistics', 'p-value']):
            return """**Statistical Significance in Microbiome Analysis:**

InfIntE ensures statistical rigor through:

1. **Permutation Testing**: Uses multiple random subsamples (nperms parameter)
2. **StARS Selection**: Stability Approach to Regularization Selection
3. **Compression Values**: Quantify interaction strength and reliability
4. **Cross-Validation**: Tests model stability across different data subsets

**Interpretation:**
- Higher compression values = stronger, more reliable interactions
- StARS threshold (default 0.01) controls false discovery rate
- Multiple permutations ensure results aren't due to chance
            """
        
        elif any(word in query_lower for word in ['step-by-step', 'analysis process', 'how does', 'workflow']):
            return self.explain_analysis_process()
        
        elif any(word in query_lower for word in ['interaction types', 'types of interactions', 'what interactions']):
            types_explanation = "**Interaction Types in InfIntE:**\n\n"
            for interaction_type, description in self.interaction_types.items():
                types_explanation += f"**{interaction_type}:** {description}\n\n"
            return types_explanation
        
        elif any(term in query_lower for term in ["status", "current analysis", "what do we have", "summary", "combined analysis"]):
            if self.current_analysis:
                interactions = self.current_analysis['interactions']
                if "summary" in query_lower or "combined" in query_lower:
                    return self.get_combined_analysis_summary()
                else:
                    return f"Current analysis contains {len(interactions)} interactions. You can ask for network visualization or specific interaction details."
            else:
                return "No analysis data available. Please upload and analyze microbiome data first."
        
        else:
            return "I'm not sure how to answer that. Type 'help' to see available commands, or ask me about microbiome interactions, analysis processes, or specific species relationships."

    def process_user_query_with_files(self, query: str) -> str:
        """Enhanced query processing that includes file content search"""
        query_lower = query.lower()
        
        # Check for PDF search queries
        if any(term in query_lower for term in ["search", "find", "paper", "document", "literature", "study"]) and self.uploaded_files['pdf_files']:
            return self.search_across_all_pdfs(query)
        
        # Fall back to regular query processing
        return self.process_user_query(query)

    def filter_interactions(self, filter_type: str = None, species: str = None, min_compression: float = None) -> Dict[str, Any]:
        """Filter interactions based on various criteria"""
        if not self.current_analysis:
            return {"error": "No analysis data available"}
        
        interactions = self.current_analysis['interactions'].copy()
        original_count = len(interactions)
        
        # Filter by interaction type
        if filter_type:
            interactions = interactions[interactions['lnk'].str.contains(filter_type, case=False, na=False)]
        
        # Filter by species
        if species:
            interactions = interactions[
                (interactions['sp1'].str.contains(species, case=False, na=False)) |
                (interactions['sp2'].str.contains(species, case=False, na=False))
            ]
        
        # Filter by minimum compression value
        if min_compression is not None:
            interactions = interactions[interactions['comp'] >= min_compression]
        
        filtered_count = len(interactions)
        
        return {
            'success': True,
            'filtered_interactions': interactions.to_dict('records'),
            'original_count': original_count,
            'filtered_count': filtered_count,
            'filters_applied': {
                'type': filter_type,
                'species': species,
                'min_compression': min_compression
            }
        }
    
    def get_interaction_summary(self) -> Dict[str, Any]:
        """Get summary statistics of current interactions"""
        if not self.current_analysis:
            return {"error": "No analysis data available"}
        
        interactions = self.current_analysis['interactions']
        
        summary = {
            'total_interactions': len(interactions),
            'interaction_types': interactions['lnk'].value_counts().to_dict(),
            'species_involved': len(set(interactions['sp1'].tolist() + interactions['sp2'].tolist())),
            'compression_stats': {
                'mean': float(interactions['comp'].mean()),
                'median': float(interactions['comp'].median()),
                'min': float(interactions['comp'].min()),
                'max': float(interactions['comp'].max())
            },
            'top_species': {
                'most_interactions': interactions['sp1'].value_counts().head(5).to_dict(),
                'most_targeted': interactions['sp2'].value_counts().head(5).to_dict()
            }
        }
        
        return summary
    
    def search_pdf_content(self, query: str) -> List[Dict[str, str]]:
        """Search through uploaded PDF content for relevant information"""
        if not self.uploaded_files['pdf_files']:
            return []
        
        query_lower = query.lower()
        results = []
        
        for filename, file_info in self.uploaded_files['pdf_files'].items():
            content = file_info['content'].lower()
            
            # Simple keyword matching - find sentences containing query terms
            sentences = content.split('.')
            relevant_sentences = []
            
            for sentence in sentences:
                sentence = sentence.strip()
                if len(sentence) > 20:  # Skip very short sentences
                    # Check if sentence contains query keywords
                    query_words = query_lower.split()
                    matches = sum(1 for word in query_words if word in sentence)
                    
                    if matches > 0:
                        relevant_sentences.append((sentence, matches))
            
            # Sort by relevance and take top matches
            relevant_sentences.sort(key=lambda x: x[1], reverse=True)
            
            if relevant_sentences:
                # Combine top sentences into context
                context_sentences = [s[0] for s in relevant_sentences[:3]]
                context = '. '.join(context_sentences)
                
                # Limit context length
                if len(context) > 500:
                    context = context[:500] + "..."
                
                results.append({
                    'filename': filename,
                    'context': context,
                    'relevance_score': sum(s[1] for s in relevant_sentences[:3])
                })
        
        # Sort results by relevance score
        results.sort(key=lambda x: x['relevance_score'], reverse=True)
        return results

    def analyze_multiple_csv_files(self, csv_files_data: List[pd.DataFrame]) -> Dict[str, Any]:
        """Analyze multiple CSV files and combine results"""
        if not csv_files_data:
            return {"error": "No CSV files provided"}
        
        combined_results = []
        total_interactions = 0
        
        for i, csv_data in enumerate(csv_files_data):
            try:
                result = self.analyze_microbiome_data(csv_data)
                if 'error' not in result:
                    combined_results.append({
                        'file_index': i + 1,
                        'interactions': result.get('num_interactions', 0),
                        'types': result.get('interaction_types', {})
                    })
                    total_interactions += result.get('num_interactions', 0)
                else:
                    combined_results.append({
                        'file_index': i + 1,
                        'error': result['error']
                    })
            except Exception as e:
                combined_results.append({
                    'file_index': i + 1,
                    'error': str(e)
                })
        
        # Generate summary message
        successful_files = [r for r in combined_results if 'error' not in r]
        failed_files = [r for r in combined_results if 'error' in r]
        
        if successful_files:
            message = f"🧬 Multi-file analysis complete! Processed {len(successful_files)} files successfully. "
            message += f"Total interactions found: {total_interactions}. "
            if failed_files:
                message += f"({len(failed_files)} files failed to process)"
        else:
            message = f"❌ All {len(failed_files)} files failed to process"
        
        return {
            'success': len(successful_files) > 0,
            'message': message,
            'total_interactions': total_interactions,
            'files_processed': len(successful_files),
            'files_failed': len(failed_files),
            'detailed_results': combined_results
        }

    def get_combined_analysis_summary(self) -> str:
        """Generate summary for combined analysis from multiple files"""
        if not self.current_analysis:
            return "No analysis data available. Please upload and analyze CSV files first."
        
        interactions = self.current_analysis['interactions']
        total_interactions = len(interactions)
        
        # Get interaction type distribution
        type_counts = interactions['lnk'].value_counts().to_dict()
        type_summary = ", ".join([f"{k}: {v}" for k, v in type_counts.items()])
        
        # Get species involvement
        all_species = set(interactions['sp1'].tolist() + interactions['sp2'].tolist())
        
        summary = f"""
**📊 Combined Analysis Summary:**

**Total Interactions:** {total_interactions}
**Species Involved:** {len(all_species)}
**Interaction Types:** {type_summary}

**Top Interacting Species:**
"""
        
        # Get most active species
        species_counts = interactions['sp1'].value_counts().head(5)
        for species, count in species_counts.items():
            summary += f"• {species}: {count} interactions\n"
        
        return summary.strip()

    def process_multiple_pdf_files(self, pdf_files_info: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Process multiple PDF files and combine their content for search"""
        if not pdf_files_info:
            return {"error": "No PDF files provided"}
        
        processed_files = []
        total_text_length = 0
        failed_files = []
        
        for pdf_info in pdf_files_info:
            try:
                filename = pdf_info['filename']
                file_path = pdf_info['file_path']
                
                # Extract text from PDF
                text = ""
                try:
                    # Try PyMuPDF first
                    doc = fitz.open(file_path)
                    for page in doc:
                        text += page.get_text()
                    doc.close()
                    
                    if not text.strip():
                        # Fallback to PyPDF2
                        with open(file_path, 'rb') as pdf_file:
                            pdf_reader = PyPDF2.PdfReader(pdf_file)
                            for page in pdf_reader.pages:
                                text += page.extract_text()
                
                except Exception as e:
                    failed_files.append({'filename': filename, 'error': str(e)})
                    continue
                
                # Store PDF info
                self.uploaded_files['pdf_files'][filename] = {
                    'path': file_path,
                    'content': text,
                    'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                }
                
                processed_files.append({
                    'filename': filename,
                    'text_length': len(text),
                    'word_count': len(text.split())
                })
                total_text_length += len(text)
                
            except Exception as e:
                failed_files.append({'filename': pdf_info.get('filename', 'unknown'), 'error': str(e)})
        
        # Generate summary message
        if processed_files:
            message = f"📄 Multi-PDF processing complete! Successfully processed {len(processed_files)} PDF files. "
            message += f"Total text extracted: {total_text_length:,} characters. "
            if failed_files:
                message += f"({len(failed_files)} files failed to process)"
            message += " You can now search across all uploaded PDFs using queries like 'search for competition in papers'."
        else:
            message = f"❌ All {len(failed_files)} PDF files failed to process"
        
        return {
            'success': len(processed_files) > 0,
            'message': message,
            'files_processed': len(processed_files),
            'files_failed': len(failed_files),
            'total_text_length': total_text_length,
            'processed_files': processed_files,
            'failed_files': failed_files
        }

    def search_across_all_pdfs(self, query: str) -> str:
        """Enhanced PDF search that works across all uploaded PDFs with better formatting"""
        if not self.uploaded_files['pdf_files']:
            return "No PDF files uploaded. Please upload PDF files first to search their content."
        
        results = self.search_pdf_content(query)
        
        if not results:
            return f"No relevant information found for '{query}' in the {len(self.uploaded_files['pdf_files'])} uploaded PDF files."
        
        response = f"🔍 **Search Results for '{query}'** (found in {len(results)} documents):\n\n"
        
        for i, result in enumerate(results[:5], 1):  # Top 5 results
            response += f"**{i}. From {result['filename']}:**\n"
            response += f"{result['context']}\n"
            response += f"*Relevance: {result['relevance_score']:.2f}*\n\n"
        
        if len(results) > 5:
            response += f"*...and {len(results) - 5} more results. Try a more specific query for better results.*"
        
        return response.strip()

# Initialize chatbot
chatbot = MicrobiomeXAIChatbot()

@app.route('/')
def index():
    """Main chatbot interface"""
    return render_template('index.html')

@app.route('/chat', methods=['POST'])
def chat():
    """Handle chat messages"""
    data = request.json
    user_message = data.get('message', '')
    
    # Process the query
    response = chatbot.process_user_query_with_files(user_message)
    
    if response == "SHUTDOWN_REQUESTED":
        shutdown_server()
        return jsonify({'response': 'Server shutting down...'})
    
    if response == "SHOW_NETWORK":
        return jsonify({
            'response': 'Generating network visualization...',
            'show_network': True
        })

    if response == "SHOW_DEPTH":
        return jsonify({
            'response': 'Generating depth plots...',
            'show_depth': True
            })
    
    # Add to chat history
    chatbot.chat_history.append({
        'user': user_message,
        'bot': response,
        'timestamp': datetime.now().isoformat()
    })
    
    return jsonify({
        'response': response,
        'timestamp': datetime.now().isoformat()
    })

@app.route('/upload_depth', methods=['POST'])
def upload_depth():
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': "No file part named 'file'"}), 400
        f = request.files['file']
        if not f.filename:
            return jsonify({'success': False, 'error': "No selected file"}), 400
        os.makedirs(os.path.dirname(DEPTH_CSV_DEFAULT), exist_ok=True)
        f.save(DEPTH_CSV_DEFAULT)
        chatbot.load_depth_csv(DEPTH_CSV_DEFAULT)
        return jsonify({'success': True, 'schema': chatbot.depth_schema})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 400

@app.route('/depth', methods=['GET'])
def depth_viewer():
    try:
        figures = chatbot.generate_depth_figures()
        return jsonify({"success": True, **figures})
    except FileNotFoundError as e:
        return jsonify({"success": False, "error": str(e)}), 404
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 400


@app.route('/analyze', methods=['POST'])
def analyze():
    """Analyze uploaded microbiome data"""
    try:
        data = request.json
        otu_data = pd.DataFrame(data['otu_data'])
        depth = data.get('depth', None)
        
        result = chatbot.analyze_microbiome_data(otu_data, depth)
        return jsonify(result)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 400

@app.route('/network')
def network():
    """Generate network visualization"""
    result = chatbot.generate_network_visualization()
    return jsonify(result)

@app.route('/history')
def history():
    """Get chat history"""
    return jsonify(chatbot.chat_history)

@app.route('/filter_interactions', methods=['POST'])
def filter_interactions():
    """Filter interactions based on user input"""
    try:
        data = request.json
        filter_type = data.get('filter_type', None)
        species = data.get('species', None)
        min_compression = data.get('min_compression', None)
        
        result = chatbot.filter_interactions(filter_type, species, min_compression)
        return jsonify(result)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 400

@app.route('/interaction_summary')
def interaction_summary():
    """Get summary statistics of current interactions"""
    result = chatbot.get_interaction_summary()
    return jsonify(result)

@app.route('/upload_file', methods=['POST'])
def upload_file():
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': 'No file provided'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No file selected'}), 400
        
        filename = secure_filename(file.filename)
        file_ext = filename.rsplit('.', 1)[1].lower() if '.' in filename else ''
        
        if file_ext not in chatbot.allowed_extensions:
            return jsonify({'success': False, 'error': f'File type .{file_ext} not allowed'}), 400
        
        # Save file
        file_path = os.path.join(chatbot.upload_dir, filename)
        file.save(file_path)
        
        # Process based on file type
        if file_ext == 'pdf':
            try:
                # Try PyMuPDF first
                doc = fitz.open(file_path)
                text = ""
                for page in doc:
                    text += page.get_text()
                doc.close()
                
                if not text.strip():
                    # Fallback to PyPDF2
                    with open(file_path, 'rb') as pdf_file:
                        pdf_reader = PyPDF2.PdfReader(pdf_file)
                        text = ""
                        for page in pdf_reader.pages:
                            text += page.extract_text()
                
                # Store PDF info
                chatbot.uploaded_files['pdf_files'][filename] = {
                    'path': file_path,
                    'content': text,
                    'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                }
                
                return jsonify({
                    'success': True,
                    'message': f'📄 PDF {filename} uploaded and processed successfully',
                    'text_length': len(text)
                })
                
            except Exception as e:
                return jsonify({'success': False, 'error': f'PDF processing failed: {str(e)}'}), 400
        
        elif file_ext == 'csv':
            try:
                df = pd.read_csv(file_path)
                chatbot.uploaded_files['csv_files'][filename] = {
                    'path': file_path,
                    'data': df,
                    'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                }
                
                # Auto-analyze CSV file
                try:
                    analysis_result = chatbot.analyze_microbiome_data(df)
                    
                    if 'error' not in analysis_result:
                        success_message = f'🧬 Analysis complete! Found {analysis_result.get("num_interactions", 0)} interactions.'
                        chatbot.chat_history.append({
                            'user': f'Uploaded {filename}',
                            'bot': success_message,
                            'timestamp': datetime.now().isoformat()
                        })
                        
                        return jsonify({
                            'success': True,
                            'message': success_message,
                            'rows': len(df),
                            'columns': len(df.columns),
                            'interactions_found': analysis_result.get('num_interactions', 0)
                        })
                    else:
                        return jsonify({
                            'success': True,
                            'message': f'📊 CSV {filename} uploaded successfully (analysis failed: {analysis_result["error"]})',
                            'rows': len(df),
                            'columns': len(df.columns)
                        })
                except Exception as analysis_error:
                    return jsonify({
                        'success': True,
                        'message': f'📊 CSV {filename} uploaded successfully (analysis failed: {str(analysis_error)})',
                        'rows': len(df),
                        'columns': len(df.columns)
                    })
                
            except Exception as e:
                return jsonify({'success': False, 'error': f'CSV processing failed: {str(e)}'}), 400
                
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/get_uploaded_files', methods=['GET'])
def get_uploaded_files():
    """Get list of uploaded files"""
    try:
        files_info = {
            'csv_files': {},
            'pdf_files': {}
        }
        
        for filename, info in chatbot.uploaded_files['csv_files'].items():
            files_info['csv_files'][filename] = {
                'timestamp': info['timestamp'],
                'rows': len(info['data']) if 'data' in info else 0
            }
        
        for filename, info in chatbot.uploaded_files['pdf_files'].items():
            files_info['pdf_files'][filename] = {
                'timestamp': info['timestamp'],
                'text_length': len(info['content']) if 'content' in info else 0
            }
        
        return jsonify({
            'success': True,
            'files': files_info
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/search_pdf_content', methods=['POST'])
def search_pdf_content():
    """Search through PDF content"""
    try:
        data = request.json
        query = data.get('query', '')
        
        if not query:
            return jsonify({'success': False, 'error': 'Query is required'}), 400
        
        results = chatbot.search_pdf_content(query)
        
        return jsonify({
            'success': True,
            'results': results,
            'query': query
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/upload_multiple_pdfs', methods=['POST'])
def upload_multiple_pdfs():
    """Handle multiple PDF file uploads"""
    try:
        if 'files' not in request.files:
            return jsonify({'success': False, 'error': 'No files provided'}), 400
        
        files = request.files.getlist('files')
        if not files or all(f.filename == '' for f in files):
            return jsonify({'success': False, 'error': 'No files selected'}), 400
        
        pdf_files_info = []
        
        for file in files:
            if file and file.filename:
                filename = secure_filename(file.filename)
                file_ext = filename.rsplit('.', 1)[1].lower() if '.' in filename else ''
                
                if file_ext != 'pdf':
                    continue  # Skip non-PDF files
                
                # Save file
                file_path = os.path.join(chatbot.upload_dir, filename)
                file.save(file_path)
                
                pdf_files_info.append({
                    'filename': filename,
                    'file_path': file_path
                })
        
        if not pdf_files_info:
            return jsonify({'success': False, 'error': 'No valid PDF files found'}), 400
        
        # Process all PDFs
        result = chatbot.process_multiple_pdf_files(pdf_files_info)
        
        # Add to chat history if successful
        if result.get('success'):
            chatbot.chat_history.append({
                'user': f'Uploaded {len(pdf_files_info)} PDF files',
                'bot': result['message'],
                'timestamp': datetime.now().isoformat()
            })
        
        return jsonify(result)
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

def signal_handler(sig, frame):
    print('\nShutting down server...')
    # Clean up R environment
    try:
        subprocess.check_output(['R', '--slave', '-e', 'rm(list=ls())'])
        subprocess.check_output(['R', '--slave', '-e', 'gc()'])
    except:
        pass
    # Force exit
    os._exit(0)

signal.signal(signal.SIGINT, signal_handler)

def shutdown_server():
    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=5000, threaded=False, use_reloader=False)
