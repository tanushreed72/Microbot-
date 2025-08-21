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
        # 1) Prefer freshly uploaded depth.csv held in memory
        depth_files = [
            f for f in self.uploaded_files.get('csv_files', {}).keys()
            if 'depth' in f.lower()
        ]
        if depth_files:
            depth_filename = depth_files[0]
            df = self.uploaded_files['csv_files'][depth_filename]['data'].copy()
            
            depth_col = None
            sample_col = None

            for col in df.columns:
                if pd.api.types.is_numeric_dtype(df[col]) and any(
                    kw in col.lower() for kw in ['depth', 'reads', 'count', 'size']
                ):
                    depth_col = col
                    break

            for col in df.columns:
                if any(kw in col.lower() for kw in ['sample', 'id', 'name']):
                    sample_col = col
                    break

            if depth_col is None:
                # Fallback: first numeric column
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) > 0:
                    depth_col = numeric_cols[0]
                else:
                    raise ValueError("No numeric depth column found in uploaded depth file")

            # Standardize column names
            df = df.rename(columns={depth_col: 'Depth'})
            if sample_col:
                df = df.rename(columns={sample_col: 'Sample'})
            else:
                if df.index.name is not None:
                    df['Sample'] = df.index.astype(str)
                else:
                    df['Sample'] = [f"S{i+1}" for i in range(len(df))]

        # 2) Else, if we already have a loaded DataFrame in memory, use it
        elif hasattr(self, 'depth_df') and isinstance(getattr(self, 'depth_df', None), pd.DataFrame) and not self.depth_df.empty:
            df = self.depth_df.copy()

        # 3) Else, try the default depth file on disk
        elif os.path.exists(DEPTH_CSV_DEFAULT):
            self.load_depth_csv(DEPTH_CSV_DEFAULT)
            df = self.depth_df.copy()

        # 4) Nothing available
        else:
            raise FileNotFoundError(
                "No depth data available. Please upload a depth.csv file or analyze microbiome data with depth information first.")

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
            
            # Create R script for analysis with optimized parameters using unique temp files
            input_path = os.path.join(os.getcwd(), 'input.csv')
            processed_data.to_csv(input_path, index=True)
            
            interactions_output = os.path.join(os.getcwd(), 'interactions.csv')
            interactions_enhanced_output = os.path.join(os.getcwd(), 'interactions_enhanced.csv')
            
            # Fix Windows compatibility: mclapply doesn't support multiple cores on Windows
            if os.name == 'nt':  # Windows
                ncores = 1
            else:
                ncores = min(6, os.cpu_count() or 1)  # Unix/Linux systems
            
            # Get optimized nperms from preprocessing if available, otherwise use default
            nperms = getattr(processed_data, '_nperms', 10)
            
            depth_str = ', '.join(map(str, depth)) if depth else ''
            
            # Convert paths to forward slashes for R compatibility
            input_path_r = input_path.replace('\\', '/')
            interactions_output_r = interactions_output.replace('\\', '/')
            interactions_enhanced_output_r = interactions_enhanced_output.replace('\\', '/')
            
            original_rows, original_cols = len(processed_data), processed_data.shape[1]
            
            # Adaptive timeout based on dataset size
            if original_rows > 5000 or original_cols > 1000:
                timeout_seconds = 1200  # 20 minutes for ultra-large datasets
                print(f"Ultra-large dataset detected: {original_rows}x{original_cols} - timeout: {timeout_seconds//60} minutes")
            elif original_rows > 1000 or original_cols > 500:
                timeout_seconds = 600   # 10 minutes for large datasets
                print(f"Large dataset detected: {original_rows}x{original_cols} - timeout: {timeout_seconds//60} minutes")
            else:
                timeout_seconds = 300   # 5 minutes for smaller datasets
                print(f"Small dataset detected: {original_rows}x{original_cols} - timeout: {timeout_seconds//60} minutes")
            
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
                Source = interactions$sp1,
                Target = interactions$sp2,
                Interaction_Type = interactions$lnk,
                Confidence = interactions$comp,
                Dataset_Size = "{original_rows}x{original_cols}",
                Optimization_Applied = "{'Yes' if original_rows > 100 or original_cols > 50 else 'No'}",
                Analysis_Date = Sys.Date()
            )
            write.csv(interactions_enhanced, "{interactions_enhanced_output_r}", row.names=FALSE)
            """
            
            # Write R script
            with open('script.R', 'w') as f:
                f.write(r_script)
            
            # Execute R script with extended timeout for large datasets
            print(f"Running R analysis with timeout of {timeout_seconds} seconds...")
            result = subprocess.run([self.r_executable, '--slave', '-f', 'script.R'], 
                                  capture_output=True, text=True, timeout=timeout_seconds, cwd=os.getcwd())  # Increased timeout
            
            if result.returncode != 0:
                return {"error": f"R analysis failed: {result.stderr}"}
            
            # Read results from unique interactions file
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
            return {"error": f"Analysis timed out after {timeout_seconds//60} minutes"}
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
            
            # Scalable optimization for very large datasets (5000+ rows/columns)
            original_rows, original_cols = len(processed_data), processed_data.shape[1]
            
            # Adaptive optimization based on dataset size
            if original_rows > 5000 or original_cols > 1000:
                # Ultra-large dataset optimization
                nperms = 2
                max_species = 100
                max_samples = 100
                print(f"Ultra-large dataset detected ({original_rows}x{original_cols})")
            elif original_rows > 1000 or original_cols > 500:
                # Large dataset optimization  
                nperms = 3
                max_species = 200
                max_samples = 150
                print(f"Large dataset detected ({original_rows}x{original_cols})")
            elif original_rows > 100 or original_cols > 50:
                # Medium dataset optimization
                nperms = 5
                max_species = 300
                max_samples = 200
                print(f"Medium dataset detected ({original_rows}x{original_cols})")
            else:
                # Small dataset - minimal optimization
                nperms = 10
                max_species = original_rows
                max_samples = original_cols
            
            # Intelligent species selection (rows)
            if len(processed_data) > max_species:
                # Use variance + abundance for better species selection
                row_sums = processed_data.sum(axis=1)
                row_vars = processed_data.var(axis=1)
                # Combined score: abundance + variance (normalized)
                abundance_norm = (row_sums - row_sums.min()) / (row_sums.max() - row_sums.min() + 1e-10)
                variance_norm = (row_vars - row_vars.min()) / (row_vars.max() - row_vars.min() + 1e-10)
                combined_score = abundance_norm + variance_norm
                top_species = combined_score.nlargest(max_species).index
                processed_data = processed_data.loc[top_species]
                print(f"Species optimization: {original_rows} → {len(processed_data)} (kept most abundant + variable)")
            
            # Intelligent sample selection (columns)
            if processed_data.shape[1] > max_samples:
                # Select samples with highest total abundance and diversity
                col_sums = processed_data.sum(axis=0)
                col_diversity = (processed_data > 0).sum(axis=0)  # Number of non-zero species
                # Combined score: abundance + diversity
                col_abundance_norm = (col_sums - col_sums.min()) / (col_sums.max() - col_sums.min() + 1e-10)
                col_diversity_norm = (col_diversity - col_diversity.min()) / (col_diversity.max() - col_diversity.min() + 1e-10)
                col_combined_score = col_abundance_norm + col_diversity_norm
                top_samples = col_combined_score.nlargest(max_samples).index
                processed_data = processed_data[top_samples]
                print(f"Sample optimization: {original_cols} → {len(processed_data.columns)} (kept most abundant + diverse)")
            
            # Memory-efficient processing for very large datasets
            if original_rows > 2000 or original_cols > 500:
                # Additional memory optimizations
                processed_data = processed_data.astype('float32')  # Reduce memory usage
                print(f"Applied memory optimization: using float32 precision")
            
            print(f"Final optimized dataset: {len(processed_data)}x{processed_data.shape[1]} with {nperms} permutations")
            
            # Store nperms value for later use
            processed_data._nperms = nperms
            
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
        """Process user query and return appropriate response"""
        query_lower = query.lower().strip()
        
        # Handle shutdown request
        if query_lower in ['shutdown', 'exit', 'quit']:
            return "SHUTDOWN_REQUESTED"
        
        # Handle static Q&A first (these should work regardless of uploaded files)
        if 'types' in query_lower and 'interaction' in query_lower:
            types_explanation = "**Interaction Types in InfIntE:**\n\n"
            for interaction_type, description in self.interaction_types.items():
                types_explanation += f"**{interaction_type}:** {description}\n\n"
            return types_explanation
        
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
- Higher compression values indicate stronger statistical support
- StARS ensures model stability across different data subsets
- Permutation testing validates statistical significance
            """
        
        elif any(word in query_lower for word in ['network topology', 'topology', 'patterns']):
            return """**Network Topology and Patterns:**

Microbiome interaction networks reveal important ecological patterns:

1. **Hub Species**: Highly connected microbes that influence many others
2. **Modules/Communities**: Groups of tightly interacting species
3. **Network Density**: Overall connectivity level in the ecosystem
4. **Centrality Measures**: Identify keystone species in the network

**Common Patterns:**
- **Scale-free networks**: Few highly connected hubs, many low-degree nodes
- **Small-world properties**: High clustering with short path lengths
- **Modular structure**: Distinct functional communities
            """
        
        elif any(word in query_lower for word in ['step-by-step', 'analysis process', 'how does', 'process']):
            return self.explain_analysis_process()
        
        # Handle depth plot requests FIRST (before network to avoid conflict with 'plot' keyword)
        # Handle depth plot requests FIRST (before network to avoid conflict)
        depth_terms = [
            "depth plot", "show depth", "sequencing depth", "depth distribution",
            "library size", "depth visualization", "depth stats", "depth qc", "show depth plot"
        ]
        if any(term in query_lower for term in depth_terms):
            return "SHOW_DEPTH"


        # Handle network visualization requests (more specific terms to avoid conflict)
        elif any(term in query_lower for term in ['network', 'show me the network', 'interactions network']) or \
             (any(term in query_lower for term in ['plot', 'graph', 'visualization']) and any(term in query_lower for term in ['interaction', 'network'])):
            if not self.current_analysis:
                return "No analysis data available. Please upload and analyze data first."
            return "SHOW_NETWORK"
        
        # Handle CSV data queries only after static Q&A
        csv_response = self.handle_csv_data_query(query_lower)
        if csv_response:
            return csv_response
        
        # ---- Current Analysis Queries (High Priority) ----
        # Dataset filtering commands - these should work BEFORE PDF search
        if any(word in query_lower for word in ['filter', 'show me', 'find', 'top', 'interactions with']) and self.current_analysis:
            try:
                if 'abundance' in query_lower or any(pattern in query_lower for pattern in ['> 0.1', '> 0.5', '>0.1', '>0.5']):
                    # Extract abundance threshold
                    import re
                    abundance_match = re.search(r'abundance\s*[><=]+\s*([0-9.]+)', query_lower)
                    threshold_match = re.search(r'>\s*([0-9.]+)', query_lower)
                    
                    if abundance_match:
                        threshold = float(abundance_match.group(1))
                    elif threshold_match:
                        threshold = float(threshold_match.group(1))
                    else:
                        threshold = 0.1  # Default threshold
                    
                    # Filter interactions by compression value (proxy for abundance)
                    interactions = self.current_analysis['interactions']
                    if 'comp' in interactions.columns:
                        # Scale threshold appropriately for compression values
                        if threshold <= 1.0:
                            scaled_threshold = threshold * 1000  # Scale 0.1 to 100, 0.5 to 500
                        else:
                            scaled_threshold = threshold
                        
                        filtered = interactions[interactions['comp'] > scaled_threshold]
                        
                        if len(filtered) > 0:
                            return f"**Found {len(filtered)} interactions with compression > {scaled_threshold:.0f} (abundance > {threshold}):**\n\n" + \
                                   filtered.nlargest(min(10, len(filtered)), 'comp')[['sp1', 'sp2', 'lnk', 'comp']].to_string(index=False)
                        else:
                            return f"No interactions found with compression > {scaled_threshold:.0f} (abundance > {threshold}). Try a lower threshold."
                    else:
                        return "Compression data not available for abundance filtering."
                            
                elif 'top' in query_lower:
                    # Handle "top N interactions" or just "top N" 
                    import re
                    num_match = re.search(r'top\s+(\d+)', query_lower)
                    if num_match:
                        n = int(num_match.group(1))
                        interactions = self.current_analysis['interactions']
                        if 'comp' in interactions.columns:
                            top_interactions = interactions.nlargest(n, 'comp')
                            
                            # Format output with better species name display
                            result = f"**Top {n} interactions by compression value:**\n\n"
                            for idx, row in top_interactions.iterrows():
                                sp1_short = str(row['sp1'])[:50] + "..." if len(str(row['sp1'])) > 50 else str(row['sp1'])
                                sp2_short = str(row['sp2'])[:50] + "..." if len(str(row['sp2'])) > 50 else str(row['sp2'])
                                result += f"**{row['lnk']}** between:\n"
                                result += f"  • {sp1_short}\n"
                                result += f"  • {sp2_short}\n"
                                result += f"  • Compression: {row['comp']}\n\n"
                            
                            return result
                        else:
                            return "Compression data not available for ranking."
                    else:
                        # Default to top 10 if no number specified
                        interactions = self.current_analysis['interactions']
                        if 'comp' in interactions.columns:
                            top_interactions = interactions.nlargest(10, 'comp')
                            
                            result = f"**Top 10 interactions by compression value:**\n\n"
                            for idx, row in top_interactions.iterrows():
                                sp1_short = str(row['sp1'])[:50] + "..." if len(str(row['sp1'])) > 50 else str(row['sp1'])
                                sp2_short = str(row['sp2'])[:50] + "..." if len(str(row['sp2'])) > 50 else str(row['sp2'])
                                result += f"**{row['lnk']}** between:\n"
                                result += f"  • {sp1_short}\n"
                                result += f"  • {sp2_short}\n"
                                result += f"  • Compression: {row['comp']}\n\n"
                            
                            return result
                        else:
                            return "Compression data not available for ranking."
            except Exception as e:
                return f"Error processing filter request: {str(e)}. Please try uploading and analyzing data first."
        
        # ---- PDF Search Queries ----
        # Check for PDF search queries ONLY with explicit PDF references
        explicit_pdf_terms = ["pdf","according to the paper", "paper", "from the pdf", "paper says", "document states", "literature shows", "study mentions", "research paper", "in pdf"]
        if any(term in query_lower for term in explicit_pdf_terms) and self.uploaded_files['pdf_files']:
            return self.search_across_all_pdfs(query)
        
        # ---- CSV Data Queries ----
        # Handle questions about uploaded CSV data
        if self.uploaded_files['csv_files']:
            csv_response = self.handle_intelligent_csv_query(query_lower)
            if csv_response:
                return csv_response
        
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
        
        # These are handled earlier in the function now
        
        elif 'help' in query_lower:
            return """
**Available Commands:**
- Ask about interaction types: "What are the interaction types?"
- Explain specific interactions: "Explain interaction between SpeciesA and SpeciesB"  
- Learn about the process: "How does the analysis work?"
- Request visualizations: "Show me the network"
- Query uploaded data: "List species", "Show taxonomy", "What phyla are present?"
- Get help: "help"

**Example Queries:**
- "What is mutualism?"
- "Explain the interaction between Bacillus and Pseudomonas"
- "How does InfIntE work?"
- "Show me interaction types"
- "List all species in the data"
- "What families are present?"
            """
        
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

    def analyze_csv_file(self, file_path: str) -> Dict[str, Any]:
        """Analyze CSV file and store in uploaded_files for query handling"""
        try:
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            # Store in uploaded_files for query handling
            filename = os.path.basename(file_path)
            self.uploaded_files['csv_files'][filename] = {
                'path': file_path,
                'data': df,
                'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            }
            
            # Detect CSV type for intelligent processing
            csv_type = self.detect_csv_type(df, filename)
            
            # Generate automatic analysis summary
            analysis_summary = self.generate_auto_analysis_summary(df, filename, csv_type, False)
            
            # Handle depth CSV files specifically
            if csv_type == 'depth':
                # Also save to the default depth location for compatibility
                os.makedirs(os.path.dirname(DEPTH_CSV_DEFAULT), exist_ok=True)
                df.to_csv(DEPTH_CSV_DEFAULT, index=False)
                
                # Load into depth_df for immediate access
                try:
                    self.load_depth_csv(DEPTH_CSV_DEFAULT)
                except Exception as e:
                    print(f"Warning: Could not load depth file into depth_df: {e}")
                
                return {
                    'success': True,
                    'message': f'✅ Depth file {filename} uploaded successfully! Ask for a depth plot to visualize sequencing depth.',
                    'file_stored': True,
                    'auto_analysis': analysis_summary,
                    'csv_type': csv_type
                }
            
            # Handle interaction CSV files directly without R analysis
            elif csv_type == 'interactions':
                # Store interaction data directly for network visualization
                interactions = df.copy()
                
                # Standardize column names for interaction data
                if 'sp1' in df.columns and 'sp2' in df.columns:
                    # Already in correct format
                    pass
                elif 'species1' in df.columns and 'species2' in df.columns:
                    # Rename columns to standard format
                    interactions = interactions.rename(columns={
                        'species1': 'sp1',
                        'species2': 'sp2',
                        'interaction': 'lnk',
                        'type': 'lnk'
                    })
                
                # Store analysis results for network visualization
                self.current_analysis = {
                    'interactions': interactions,
                    'otu_data': df,
                    'timestamp': datetime.now().isoformat(),
                    'source': 'interaction_csv'
                }
                
                return {
                    'success': True,
                    'message': f'✅ Interaction data loaded successfully from {filename}',
                    'interactions_count': len(interactions),
                    'file_stored': True,
                    'auto_analysis': analysis_summary,
                    'csv_type': csv_type,
                    'show_network': True  # Trigger network visualization
                }
            
            # Check if this is microbiome data that needs InfIntE analysis
            is_microbiome = self.detect_microbiome_data(df)
            
            if is_microbiome:
                # Run microbiome analysis automatically
                try:
                    result = self.analyze_microbiome_data(df)
                    result['file_stored'] = True
                    result['auto_analysis'] = analysis_summary
                    result['csv_type'] = csv_type
                    
                    # Store analysis results for future queries
                    self.uploaded_files['csv_files'][filename]['analysis_results'] = result
                    
                    return result
                except Exception as e:
                    # If microbiome analysis fails, fall back to generic analysis
                    return self.run_generic_csv_analysis(df, filename, csv_type, analysis_summary)
            else:
                # For non-microbiome CSV, run generic analysis
                return self.run_generic_csv_analysis(df, filename, csv_type, analysis_summary)
                
        except Exception as e:
            return {
                'success': False,
                'error': f'Failed to analyze CSV file: {str(e)}'
            }

    def generate_auto_analysis_summary(self, df: pd.DataFrame, filename: str, csv_type: str, is_microbiome: bool) -> str:
        """Generate automatic analysis summary for uploaded CSV"""
        summary_parts = []
        
        # Basic file info
        summary_parts.append(f"📊 **{filename}** ({csv_type.title()} Data)")
        summary_parts.append(f"• **Dimensions:** {len(df)} rows × {len(df.columns)} columns")
        
        # Type-specific insights
        if csv_type == 'taxonomy':
            if 'species' in df.columns:
                unique_species = df['species'].nunique()
                summary_parts.append(f"• **Species:** {unique_species} unique species identified")
            if 'phylum' in df.columns:
                unique_phyla = df['phylum'].nunique()
                summary_parts.append(f"• **Phyla:** {unique_phyla} different phyla")
                
        elif csv_type == 'interactions':
            if 'sp1' in df.columns and 'sp2' in df.columns:
                unique_species = set(df['sp1'].unique()) | set(df['sp2'].unique())
                summary_parts.append(f"• **Interactions:** {len(df)} interactions between {len(unique_species)} species")
            if 'lnk' in df.columns:
                interaction_types = df['lnk'].value_counts()
                summary_parts.append(f"• **Types:** {', '.join([f'{k}({v})' for k, v in interaction_types.head(3).items()])}")
                
        elif csv_type == 'otu_table':
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            if len(numeric_cols) > 1:
                summary_parts.append(f"• **Samples:** {len(numeric_cols)} samples")
                summary_parts.append(f"• **Features:** {len(df)} OTUs/ASVs")
                
        elif csv_type == 'species_mapping':
            if 'species_name' in df.columns:
                unique_species = df['species_name'].nunique()
                summary_parts.append(f"• **Mapped Species:** {unique_species} species with identifiers")
        
        # Analysis capability
        if is_microbiome:
            summary_parts.append("🧬 **Microbiome Analysis:** Running InfIntE interaction inference...")
        else:
            summary_parts.append("📈 **Generic Analysis:** Correlation-based analysis available")
            
        return "\n".join(summary_parts)

    def run_generic_csv_analysis(self, df: pd.DataFrame, filename: str, csv_type: str, analysis_summary: str) -> Dict[str, Any]:
        """Run generic analysis for non-microbiome CSV files"""
        try:
            # Basic statistics
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            stats = {}
            
            if len(numeric_cols) > 1:
                # Calculate correlations for numeric data
                corr_matrix = df[numeric_cols].corr()
                strong_correlations = []
                
                for i in range(len(corr_matrix.columns)):
                    for j in range(i+1, len(corr_matrix.columns)):
                        corr_val = corr_matrix.iloc[i, j]
                        if abs(corr_val) > 0.3:  # Significant correlation
                            strong_correlations.append({
                                'var1': corr_matrix.columns[i],
                                'var2': corr_matrix.columns[j],
                                'correlation': corr_val,
                                'strength': 'Strong' if abs(corr_val) > 0.7 else 'Moderate'
                            })
                
                stats['correlations'] = strong_correlations[:20]  # Top 20
                stats['numeric_summary'] = df[numeric_cols].describe().to_dict()
            
            # Store analysis results
            self.uploaded_files['csv_files'][filename]['analysis_results'] = {
                'csv_type': csv_type,
                'statistics': stats,
                'summary': analysis_summary
            }
            
            return {
                'success': True,
                'message': f'✅ {filename} analyzed successfully',
                'auto_analysis': analysis_summary,
                'csv_type': csv_type,
                'statistics': stats,
                'file_stored': True,
                'rows': len(df),
                'columns': len(df.columns)
            }
            
        except Exception as e:
            return {
                'success': False,
                'error': f'Generic analysis failed: {str(e)}'
            }
    
    def detect_csv_type(self, df: pd.DataFrame, filename: str) -> str:
        """Detect the type of CSV file based on columns and filename"""
        filename_lower = filename.lower()
        columns_lower = [col.lower() for col in df.columns]
        
        # Check filename patterns first
        if 'taxonomy' in filename_lower:
            return 'taxonomy'
        elif 'interaction' in filename_lower:
            return 'interactions'
        elif 'otu' in filename_lower or 'asv' in filename_lower:
            return 'otu_table'
        elif 'species' in filename_lower and 'mapping' in filename_lower:
            return 'species_mapping'
        elif 'depth' in filename_lower:
            return 'depth'
        
        # Check column patterns
        if any(col in columns_lower for col in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
            return 'taxonomy'
        elif any(col in columns_lower for col in ['sp1', 'sp2', 'lnk', 'comp']) or \
             (any(col in columns_lower for col in ['species1', 'species2', 'interaction', 'type'])):
            return 'interactions'
        elif len(df.columns) > 10 and df.select_dtypes(include=[np.number]).shape[1] > 5:
            return 'otu_table'  # Likely OTU table with many numeric columns
        
        return 'generic'
    
    def detect_microbiome_data(self, df: pd.DataFrame) -> bool:
        """Detect if CSV contains microbiome data that needs InfIntE analysis"""
        columns_lower = [col.lower() for col in df.columns]
        
        # Check for OTU/ASV patterns (abundance data)
        otu_patterns = ['otu', 'asv', 'abundance', 'count']
        taxonomy_patterns = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        
        # If it has taxonomy columns, it's likely microbiome data
        if any(pattern in ' '.join(columns_lower) for pattern in taxonomy_patterns):
            return True
            
        # If it has OTU/ASV patterns and many numeric columns, it's likely abundance data
        if any(pattern in ' '.join(columns_lower) for pattern in otu_patterns):
            numeric_cols = df.select_dtypes(include=[np.number]).shape[1]
            if numeric_cols > 3:  # More than 3 numeric columns suggests abundance matrix
                return True
        
        # If it has many numeric columns (>10) and few text columns, likely OTU table
        if len(df.columns) > 10 and df.select_dtypes(include=[np.number]).shape[1] > 8:
            return True
            
        return False
    
    def handle_csv_data_query(self, query_lower: str) -> str:
        """Handle queries about uploaded CSV data dynamically - works with all CSV types"""
        try:
            # Get the most recent CSV file data
            if not self.uploaded_files['csv_files']:
                return None
            
            # Get the most recently uploaded CSV file
            latest_file = max(self.uploaded_files['csv_files'].items(), 
                            key=lambda x: x[1]['timestamp'])
            filename, file_info = latest_file
            df = file_info['data']
            
            # Detect CSV type and handle accordingly
            csv_type = self.detect_csv_type(df, filename)
            
            # Species/Taxa queries - works for all CSV types
            if any(term in query_lower for term in ['list species', 'show species', 'what species', 'species list']):
                return self.handle_species_query(df, filename, csv_type)

            # Unique names queries - for species mapping files
            elif any(term in query_lower for term in ['unique names', 'unique name', 'with its unique names', 'with unique names']):
                return self.handle_species_query(df, filename, csv_type, query_type='unique_names')
            
            # Taxonomy-specific queries with intelligent filtering
            elif csv_type == 'taxonomy':
                return self.handle_intelligent_taxonomy_query(df, filename, query_lower)
            
            # Interactions-specific queries
            elif csv_type == 'interactions':
                if any(term in query_lower for term in ['interaction types', 'types of interactions', 'what interactions']):
                    return self.handle_interaction_types_query(df, filename)
                elif any(term in query_lower for term in ['strongest interactions', 'top interactions', 'best interactions']):
                    return self.handle_top_interactions_query(df, filename)
                elif any(term in query_lower for term in ['species involved', 'which species', 'organisms']):
                    return self.handle_species_in_interactions_query(df, filename)
            
            # OTU table queries
            elif csv_type == 'otu_table':
                if any(term in query_lower for term in ['samples', 'sample names', 'what samples']):
                    return self.handle_samples_query(df, filename)
                elif any(term in query_lower for term in ['otus', 'otu list', 'what otus']):
                    return self.handle_otus_query(df, filename)
                elif any(term in query_lower for term in ['abundance', 'most abundant', 'highest abundance']):
                    return self.handle_abundance_query(df, filename)
            
            # General data queries - works for all CSV types
            elif any(term in query_lower for term in ['columns', 'what columns', 'column names']):
                return f"**Columns in {filename}:**\n\n" + \
                       "\n".join([f"• **{col}** ({df[col].dtype})" for col in df.columns])
            
            elif any(term in query_lower for term in ['rows', 'how many rows', 'row count', 'size']):
                return f"**Dataset Info for {filename}:**\n\n" + \
                       f"• **Rows**: {len(df)}\n" + \
                       f"• **Columns**: {len(df.columns)}\n" + \
                       f"• **CSV Type**: {csv_type.title()}\n" + \
                       f"• **Memory Usage**: {df.memory_usage(deep=True).sum() / 1024:.1f} KB"
            
            elif any(term in query_lower for term in ['summary', 'describe', 'overview', 'info']):
                return self.handle_data_summary_query(df, filename, csv_type)
            
            # Search functionality
            elif any(term in query_lower for term in ['find', 'search']):
                # Extract search term
                search_terms = query_lower.replace('find', '').replace('search', '').strip()
                if search_terms:
                    return self.handle_search_query(df, filename, search_terms)
                else:
                    return "Please specify what you want to search for. Example: 'find Bacillus'"
            
            return None  # No match found
            
        except Exception as e:
            return f"Error processing CSV query: {str(e)}"
    
    def handle_intelligent_taxonomy_query(self, df: pd.DataFrame, filename: str, query_lower: str) -> str:
        """Handle taxonomy queries with intelligent parsing and filtering"""
        import re
        
        # Available taxonomy columns (case-insensitive mapping)
        taxonomy_levels = {}
        for col in df.columns:
            col_lower = col.lower()
            if col_lower in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                taxonomy_levels[col_lower] = col
        
        # Extract specific taxonomic names from query
        def extract_taxonomic_names(query):
            # Look for capitalized words that might be taxonomic names
            taxonomic_names = re.findall(r'\b([A-Z][a-z]+(?:[a-z]*[A-Z]*[a-z]*)*)\b', query)
            return [name for name in taxonomic_names if len(name) > 2]
        
        # Parse relationship queries like "order X has which family"
        relationship_patterns = [
            (r'(\w+)\s+(\w+)\s+has\s+which\s+(\w+)', 'filter_by_parent'),  # "order Micrococcales has which family"
            (r'which\s+(\w+)\s+are\s+in\s+(\w+)\s+(\w+)', 'filter_by_parent'),  # "which families are in order Micrococcales"
            (r'(\w+)\s+in\s+(\w+)\s+(\w+)', 'filter_by_parent'),  # "families in order Micrococcales"
            (r'show\s+(\w+)\s+for\s+(\w+)\s+(\w+)', 'filter_by_parent'),  # "show families for order Micrococcales"
        ]
        
        for pattern, query_type in relationship_patterns:
            match = re.search(pattern, query_lower)
            if match:
                if query_type == 'filter_by_parent':
                    groups = match.groups()
                    if len(groups) == 3:
                        # Pattern: "order Micrococcales has which family"
                        parent_level = groups[0].lower()  # order
                        parent_name = groups[1]  # Micrococcales
                        child_level = groups[2].lower()  # family
                        
                        # Check if we have the required columns
                        if parent_level in taxonomy_levels and child_level in taxonomy_levels:
                            parent_col = taxonomy_levels[parent_level]
                            child_col = taxonomy_levels[child_level]
                            
                            # Filter data by parent taxonomic level
                            filtered_df = df[df[parent_col].str.contains(parent_name, case=False, na=False)]
                            
                            if len(filtered_df) > 0:
                                child_taxa = filtered_df[child_col].dropna().unique()
                                child_taxa = [t for t in child_taxa if str(t) != 'nan' and str(t).strip()]
                                
                                if len(child_taxa) > 0:
                                    return f"**{child_level.title()}s in {parent_level} '{parent_name}':**\n\n" + \
                                           "\n".join([f"• {taxon}" for taxon in sorted(child_taxa)])
                                else:
                                    return f"No {child_level}s found for {parent_level} '{parent_name}' in {filename}"
                            else:
                                return f"No records found for {parent_level} '{parent_name}' in {filename}"
        
        # Handle reverse queries like "which order has family Microbacteriaceae"
        reverse_patterns = [
            (r'which\s+(\w+)\s+has\s+(\w+)\s+(\w+)', 'find_parent'),  # "which order has family Microbacteriaceae"
            (r'(\w+)\s+(\w+)\s+belongs\s+to\s+which\s+(\w+)', 'find_parent'),  # "family Microbacteriaceae belongs to which order"
        ]
        
        for pattern, query_type in reverse_patterns:
            match = re.search(pattern, query_lower)
            if match:
                if query_type == 'find_parent':
                    groups = match.groups()
                    if len(groups) == 3:
                        parent_level = groups[0].lower()  # order
                        child_level = groups[1].lower()  # family
                        child_name = groups[2]  # Microbacteriaceae
                        
                        if parent_level in taxonomy_levels and child_level in taxonomy_levels:
                            parent_col = taxonomy_levels[parent_level]
                            child_col = taxonomy_levels[child_level]
                            
                            # Find parent for specific child
                            filtered_df = df[df[child_col].str.contains(child_name, case=False, na=False)]
                            
                            if len(filtered_df) > 0:
                                parent_taxa = filtered_df[parent_col].dropna().unique()
                                parent_taxa = [t for t in parent_taxa if str(t) != 'nan' and str(t).strip()]
                                
                                if len(parent_taxa) > 0:
                                    return f"**{parent_level.title()}s containing {child_level} '{child_name}':**\n\n" + \
                                           "\n".join([f"• {taxon}" for taxon in sorted(parent_taxa)])
                                else:
                                    return f"No {parent_level}s found containing {child_level} '{child_name}' in {filename}"
                            else:
                                return f"No records found for {child_level} '{child_name}' in {filename}"
        
        # Handle simple listing queries
        simple_queries = {
            'phyla': 'phylum', 'phylum': 'phylum', 'what phyla': 'phylum',
            'families': 'family', 'family': 'family', 'what families': 'family', 'show families': 'family',
            'classes': 'class', 'class': 'class', 'what classes': 'class',
            'orders': 'order', 'order': 'order', 'what orders': 'order',
            'genus': 'genus', 'genera': 'genus', 'what genus': 'genus', 'list genus': 'genus', 'show genus': 'genus',
            'species': 'species', 'what species': 'species', 'list species': 'species'
        }
        
        for query_term, level in simple_queries.items():
            if query_term in query_lower:
                if level in taxonomy_levels:
                    return self.handle_taxonomy_query(df, filename, level)
        
        # Handle taxonomy summary
        if any(term in query_lower for term in ['taxonomy', 'taxonomic', 'classification']):
            return self.handle_taxonomy_summary(df, filename)
        
        # Default fallback
        return f"I couldn't understand your taxonomy query. Try asking:\n" + \
               "• 'order Micrococcales has which family'\n" + \
               "• 'which order has family Microbacteriaceae'\n" + \
               "• 'show families' or 'list genus'\n" + \
               f"Available taxonomy levels: {', '.join(taxonomy_levels.keys())}"
    
    def handle_species_query(self, df: pd.DataFrame, filename: str, csv_type: str, query_type: str = None) -> str:
        """Handle species listing queries for different CSV types"""
        if csv_type == 'taxonomy':
            # Look for genus/species columns
            for col in ['genus', 'species', 'Genus', 'Species']:
                if col in df.columns:
                    species_list = df[col].dropna().unique()
                    species_list = [s for s in species_list if str(s) != 'nan' and str(s).strip()]
                    if len(species_list) > 50:
                        return f"**Species/Genera in {filename}:** (showing first 50 of {len(species_list)})\n\n" + \
                               "\n".join([f"• {species}" for species in sorted(species_list)[:50]]) + \
                               f"\n\n*...and {len(species_list) - 50} more.*"
                    else:
                        return f"**Species/Genera in {filename}:**\n\n" + \
                               "\n".join([f"• {species}" for species in sorted(species_list)])
        
        elif csv_type == 'interactions':
            # Get unique species from sp1 and sp2 columns
            species_set = set()
            for col in ['sp1', 'sp2', 'species1', 'species2']:
                if col in df.columns:
                    species_set.update(df[col].dropna().unique())
            
            if species_set:
                species_list = [s for s in species_set if str(s) != 'nan' and str(s).strip()]
                return f"**Species in Interactions ({filename}):**\n\n" + \
                       "\n".join([f"• {species}" for species in sorted(species_list)])
        
        elif csv_type == 'otu_table':
            # For OTU tables, show sample names or OTU IDs
            if 'sample' in df.columns[0].lower() or df.index.name and 'sample' in df.index.name.lower():
                samples = df.iloc[:, 0].unique() if 'sample' in df.columns[0].lower() else df.index.unique()
                return f"**Samples in {filename}:**\n\n" + \
                       "\n".join([f"• {sample}" for sample in sorted(map(str, samples))[:50]])
            else:
                # Show OTU/column names
                otus = [col for col in df.columns if not col.lower() in ['sample', 'group']][:50]
                return f"**OTUs/Features in {filename}:** (showing first 50)\n\n" + \
                       "\n".join([f"• {otu}" for otu in otus])
        
        # Generic approach for other CSV types
        potential_species_cols = [col for col in df.columns if any(term in col.lower() 
                                for term in ['species', 'genus', 'organism', 'taxa', 'otu', 'name'])]
        
        # Special handling for species mapping files
        if query_type == 'unique_names' and 'unique_name' in df.columns:
            unique_names = df['unique_name'].dropna().unique()
            unique_names = [s for s in unique_names if str(s) != 'nan' and str(s).strip()]
            if len(unique_names) > 50:
                return f"**Unique Species Names in {filename}:** (showing first 50 of {len(unique_names)})\n\n" + \
                       "\n".join([f"• {name}" for name in sorted(unique_names)[:50]]) + \
                       f"\n\n*...and {len(unique_names) - 50} more.*"
            else:
                return f"**Unique Species Names in {filename}:**\n\n" + \
                       "\n".join([f"• {name}" for name in sorted(unique_names)])
        
        elif 'species_name' in df.columns:
            species_list = df['species_name'].dropna().unique()
            species_list = [s for s in species_list if str(s) != 'nan' and str(s).strip()]
            if len(species_list) > 50:
                return f"**Species in {filename}:** (showing first 50 of {len(species_list)})\n\n" + \
                       "\n".join([f"• {species}" for species in sorted(species_list)[:50]]) + \
                       f"\n\n*...and {len(species_list) - 50} more.*"
            else:
                return f"**Species in {filename}:**\n\n" + \
                       "\n".join([f"• {species}" for species in sorted(species_list)])
        
        elif 'unique_name' in df.columns:
            unique_names = df['unique_name'].dropna().unique()
            unique_names = [s for s in unique_names if str(s) != 'nan' and str(s).strip()]
            if len(unique_names) > 50:
                return f"**Unique Species Names in {filename}:** (showing first 50 of {len(unique_names)})\n\n" + \
                       "\n".join([f"• {name}" for name in sorted(unique_names)[:50]]) + \
                       f"\n\n*...and {len(unique_names) - 50} more.*"
            else:
                return f"**Unique Species Names in {filename}:**\n\n" + \
                       "\n".join([f"• {name}" for name in sorted(unique_names)])
        
        return f"No species/taxa column clearly identified in {filename}. Available columns: {', '.join(df.columns)}"
    
    def handle_taxonomy_query(self, df: pd.DataFrame, filename: str, level: str) -> str:
        """Handle taxonomy level queries"""
        level_col = None
        for col in df.columns:
            col_lower = col.lower()
            if col_lower == level.lower():
                level_col = col
                break
        
        if level_col:
            taxa = df[level_col].dropna().unique()
            taxa = [t for t in taxa if str(t) != 'nan' and str(t).strip()]
            if len(taxa) > 30:
                return f"**{level.title()}s in {filename}:** (showing first 30 of {len(taxa)})\n\n" + \
                       "\n".join([f"• {taxon}" for taxon in sorted(taxa)[:30]]) + \
                       f"\n\n*...and {len(taxa) - 30} more.*"
            else:
                return f"**{level.title()}s in {filename}:**\n\n" + \
                       "\n".join([f"• {taxon}" for taxon in sorted(taxa)])
        else:
            return f"No {level} column found in {filename}. Available columns: {', '.join(df.columns)}"
    
    def handle_taxonomy_summary(self, df: pd.DataFrame, filename: str) -> str:
        """Handle taxonomy summary queries"""
        summary = f"**Taxonomic Summary for {filename}:**\n\n"
        for col in df.columns:
            col_lower = col.lower()
            if col_lower in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                unique_count = df[col].dropna().nunique()
                summary += f"• **{col}**: {unique_count} unique entries\n"
        return summary
    
    def handle_interaction_types_query(self, df: pd.DataFrame, filename: str) -> str:
        """Handle interaction types queries"""
        if 'lnk' in df.columns:
            types = df['lnk'].value_counts()
            return f"**Interaction Types in {filename}:**\n\n" + \
                   "\n".join([f"• **{itype}**: {count} interactions" for itype, count in types.items()])
        elif 'type' in df.columns:
            types = df['type'].value_counts()
            return f"**Interaction Types in {filename}:**\n\n" + \
                   "\n".join([f"• **{itype}**: {count} interactions" for itype, count in types.items()])
        else:
            return f"No interaction type column found in {filename}. Available columns: {', '.join(df.columns)}"
    
    def handle_top_interactions_query(self, df: pd.DataFrame, filename: str) -> str:
        """Handle top interactions queries"""
        if 'comp' in df.columns:
            top_interactions = df.nlargest(10, 'comp')[['sp1', 'sp2', 'lnk', 'comp']]
            return f"**Top 10 Interactions by Compression in {filename}:**\n\n" + \
                   top_interactions.to_string(index=False)
        else:
            return f"No compression column found in {filename} to rank interactions."
    
    def handle_species_in_interactions_query(self, df: pd.DataFrame, filename: str) -> str:
        """Handle species involved in interactions queries"""
        species_count = {}
        for col in ['sp1', 'sp2']:
            if col in df.columns:
                for species in df[col].dropna():
                    species_count[species] = species_count.get(species, 0) + 1
        
        if species_count:
            sorted_species = sorted(species_count.items(), key=lambda x: x[1], reverse=True)
            return f"**Species Involvement in {filename}:**\n\n" + \
                   "\n".join([f"• **{species}**: {count} interactions" for species, count in sorted_species[:20]])
        else:
            return f"No species columns found in {filename}."
    
    def handle_samples_query(self, df: pd.DataFrame, filename: str) -> str:
        """Handle sample queries for OTU tables"""
        if df.index.name and 'sample' in df.index.name.lower():
            samples = df.index.unique()
        elif 'sample' in df.columns[0].lower():
            samples = df.iloc[:, 0].unique()
        else:
            samples = df.index.unique()
        
        return f"**Samples in {filename}:**\n\n" + \
               "\n".join([f"• {sample}" for sample in sorted(map(str, samples))[:50]])
    
    def handle_otus_query(self, df: pd.DataFrame, filename: str) -> str:
        """Handle OTU queries"""
        otus = [col for col in df.columns if not col.lower() in ['sample', 'group']]
        return f"**OTUs/Features in {filename}:** (showing first 50 of {len(otus)})\n\n" + \
               "\n".join([f"• {otu}" for otu in otus[:50]])
    
    def handle_abundance_query(self, df: pd.DataFrame, filename: str) -> str:
        """Handle abundance queries"""
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            # Find most abundant OTUs/features
            col_sums = df[numeric_cols].sum().sort_values(ascending=False)
            return f"**Most Abundant Features in {filename}:**\n\n" + \
                   "\n".join([f"• **{col}**: {val:.2f}" for col, val in col_sums.head(20).items()])
        else:
            return f"No numeric abundance data found in {filename}."
    
    def handle_data_summary_query(self, df: pd.DataFrame, filename: str, csv_type: str) -> str:
        """Handle data summary queries"""
        summary = f"**Data Summary for {filename}:**\n\n"
        summary += f"• **File Type**: {csv_type.title()}\n"
        summary += f"• **Rows**: {len(df)}\n"
        summary += f"• **Columns**: {len(df.columns)}\n"
        summary += f"• **Numeric Columns**: {len(df.select_dtypes(include=[np.number]).columns)}\n"
        summary += f"• **Text Columns**: {len(df.select_dtypes(include=['object']).columns)}\n"
        summary += f"• **Missing Values**: {df.isnull().sum().sum()}\n"
        summary += f"• **Memory Usage**: {df.memory_usage(deep=True).sum() / 1024:.1f} KB\n\n"
        
        summary += "**Column Details:**\n"
        for col in df.columns[:10]:  # Show first 10 columns
            dtype = df[col].dtype
            non_null = df[col].count()
            summary += f"• **{col}** ({dtype}): {non_null} non-null values\n"
        
        if len(df.columns) > 10:
            summary += f"• *...and {len(df.columns) - 10} more columns*\n"
        
        return summary
    
    def handle_search_query(self, df: pd.DataFrame, filename: str, search_terms: str) -> str:
        """Handle search queries within CSV data"""
        results = []
        search_terms = search_terms.strip()
        
        # Search in text columns
        text_cols = df.select_dtypes(include=['object']).columns
        for col in text_cols:
            matches = df[df[col].astype(str).str.contains(search_terms, case=False, na=False)]
            if not matches.empty:
                results.append(f"**{col}:** {len(matches)} matches")
                # Show first few matches
                for idx, row in matches.head(5).iterrows():
                    results.append(f"  • {row[col]}")
                if len(matches) > 5:
                    results.append(f"  • *...and {len(matches) - 5} more matches*")
        
        if results:
            return f"**Search Results for '{search_terms}' in {filename}:**\n\n" + "\n".join(results)
        else:
            return f"No matches found for '{search_terms}' in {filename}."
    
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
        type_counts = interactions['lnk'].value_counts()
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
        species_counts = interactions['sp1'].value_counts()
        for species, count in species_counts.head(5).items():
            summary += f"• {species}: {count} interactions\n"
        
        return summary.strip()

    def process_pdf_file(self, filename: str, file_path: str) -> Dict[str, Any]:
        """Process a single PDF file and extract its content for search"""
        try:
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
                return {'success': False, 'error': f'Failed to extract text from PDF: {str(e)}'}
        
            if not text.strip():
                return {'success': False, 'error': 'No text could be extracted from the PDF'}
        
        # Store PDF info
            self.uploaded_files['pdf_files'][filename] = {
                'path': file_path,
                'content': text,
                'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            }
        
        # Generate success message
            word_count = len(text.split())
            message = f"📄 PDF processing complete! Successfully processed '{filename}'. "
            message += f"Extracted {len(text):,} characters ({word_count:,} words). "
            message += "You can now search this PDF using queries like 'search for microbial interactions'."
        
            return {
                'success': True,
                'message': message,
                'text_length': len(text),
                'word_count': word_count,
                'filename': filename
            }
        
        except Exception as e:
            return {'success': False, 'error': f'Error processing PDF file: {str(e)}'}

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

    def handle_intelligent_csv_query(self, query_lower: str) -> str:
        """Intelligently handle queries about uploaded CSV files"""
        try:
            # Check for OTU ID lookup queries first
            otu_lookup_result = self.handle_otu_id_lookup(query_lower)
            if otu_lookup_result:
                return otu_lookup_result
            
            # Get the most relevant CSV file for the query
            relevant_file = self.find_most_relevant_csv(query_lower)
            if not relevant_file:
                return "No CSV files available for analysis. Please upload CSV data first."
            
            filename, file_info = relevant_file
            df = file_info['data']
            csv_type = file_info.get('analysis_results', {}).get('csv_type', 'generic')
            
            # Try to answer any question about the CSV data intelligently
            return self.handle_intelligent_data_query(df, filename, csv_type, query_lower)
            
        except Exception as e:
            return f"Error processing query: {str(e)}"
    
    def handle_otu_id_lookup(self, query_lower: str) -> str:
        """Handle OTU ID to species name lookup queries and species to OTU ID queries"""
        import re

        # Check for species name to OTU ID queries FIRST (before numeric ID patterns)
        species_patterns = [
            r'otu.{0,5}id.{0,5}for.{0,5}(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'what.{0,5}is.{0,5}the.{0,5}otu.{0,5}id.{0,5}for.{0,5}(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'find.{0,5}otu.{0,5}id.{0,5}(?:for.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)'
        ]
        
        for pattern in species_patterns:
            match = re.search(pattern, query_lower)
            if match:
                species_name = match.group(1).strip()
                return self.lookup_species_to_otu(species_name)
        
        
        # Check for numeric ID queries first
        numeric_id_patterns = [
            r'(?:numeric.{0,5}id|id).{0,10}(?:of.{0,5}|for.{0,5})?(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'what.{0,5}is.{0,5}the.{0,5}(?:numeric.{0,5})?id.{0,5}(?:of.{0,5}|for.{0,5})?(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'show.{0,5}(?:numeric.{0,5})?id.{0,5}(?:of.{0,5}|for.{0,5})?(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'find.{0,5}(?:numeric.{0,5})?id.{0,5}(?:of.{0,5}|for.{0,5})?(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)'
        ]
        
        for pattern in numeric_id_patterns:
            match = re.search(pattern, query_lower)
            if match:
                species_name = match.group(1).strip()
                return self.lookup_numeric_id_for_species(species_name)
        
        # Check for OTU sequence queries first (before taxonomy patterns)
        otu_sequence_patterns = [
            r'species\s+([ACGT]{50,})\s+has\s+which\s+(order|family|class|phylum|kingdom)',
            r'([ACGT]{50,})\s+(?:has\s+which|what\s+is\s+the)\s+(order|family|class|phylum|kingdom)',
            r'what\s+(?:order|family|class|phylum|kingdom)\s+(?:does\s+)?(?:species\s+)?([ACGT]{50,})\s+(?:have|belong)'
        ]
        
        for pattern in otu_sequence_patterns:
            match = re.search(pattern, query_lower)
            if match:
                if len(match.groups()) == 2:
                    otu_sequence = match.group(1).strip().upper()
                    rank_requested = match.group(2).strip().lower()
                    return self.lookup_taxonomy_by_otu_sequence(otu_sequence, rank_requested)
        
        # Check for taxonomy queries (genus to family, family to genus, etc.)
        taxonomy_patterns = [
            r'(?:what\s+is\s+the\s+)?family\s+name\s+of\s+genus\s+([A-Za-z][A-Za-z0-9_\-]+)',
            r'(?:what\s+is\s+the\s+)?genus\s+of\s+(?:the\s+)?family\s+(?:name\s+)?([A-Za-z][A-Za-z0-9_\-]+)',
            r'(?:what\s+is\s+the\s+)?phylum\s+of\s+(?:genus\s+|order\s+)?([A-Za-z][A-Za-z0-9_\-]+)',
            r'(?:what\s+is\s+the\s+)?class\s+of\s+(?:genus\s+|order\s+)?([A-Za-z][A-Za-z0-9_\-]+)',
            r'(?:what\s+is\s+the\s+)?order\s+of\s+(?:genus\s+|phylum\s+)?([A-Za-z][A-Za-z0-9_\-]+)',
            r'(?:what\s+is\s+the\s+)?kingdom\s+of\s+(?:genus\s+|species\s+)?([A-Za-z][A-Za-z0-9_\-]+)',
            r'(?:what\s+is\s+the\s+)?species\s+of\s+(?:genus\s+)?([A-Za-z][A-Za-z0-9_\-]+)',
            r'(?:which\s+|what\s+)?species\s+(?:has\s+|have\s+|in\s+)?(?:order\s+|family\s+|class\s+|phylum\s+|kingdom\s+)([A-Za-z][A-Za-z0-9_\-]+)',
            r'(?:list\s+|show\s+|find\s+)?species\s+(?:in\s+|with\s+|from\s+)?(?:order\s+|family\s+|class\s+|phylum\s+|kingdom\s+)([A-Za-z][A-Za-z0-9_\-]+)'
        ]
        
        for pattern in taxonomy_patterns:
            match = re.search(pattern, query_lower)
            if match:
                search_term = match.group(1).strip()
                return self.lookup_taxonomy_info(search_term, query_lower)
        
        # Check for species name to OTU ID queries
        species_patterns = [
            r'otu.{0,5}id.{0,5}for.{0,5}(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'what.{0,5}is.{0,5}the.{0,5}otu.{0,5}id.{0,5}for.{0,5}(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'find.{0,5}otu.{0,5}id.{0,5}(?:for.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)'
        ]
        
        # Check for unique_name queries for specific species
        unique_name_patterns = [
            r'unique.{0,5}name.{0,5}(?:of.{0,5})?(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'what.{0,5}is.{0,5}the.{0,5}unique.{0,5}name.{0,5}(?:of.{0,5})?(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)'
        ]
        
        for pattern in unique_name_patterns:
            match = re.search(pattern, query_lower)
            if match:
                species_name = match.group(1).strip()
                return self.lookup_unique_name_for_species(species_name)
        
        for pattern in species_patterns:
            match = re.search(pattern, query_lower)
            if match:
                species_name = match.group(1).strip()
                return self.lookup_species_to_otu(species_name)
        
        # Look for long DNA sequences that might be OTU IDs
        otu_pattern = r'[ACGT]{50,}'  # DNA sequences 50+ characters long
        otu_matches = re.findall(otu_pattern, query_lower.upper())
        
        if not otu_matches:
            return None
        
        otu_id = otu_matches[0]  # Take the first/longest match
        
        # Search through uploaded CSV files for species mapping
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            # Check if this looks like a species mapping file
            if 'otu_id' in df.columns and 'species_name' in df.columns:
                # Look for exact match
                matching_rows = df[df['otu_id'].str.upper() == otu_id]
                
                if len(matching_rows) > 0:
                    species_name = matching_rows.iloc[0]['species_name']
                    return f"**OTU ID Lookup Result:**\n\n" + \
                           f"🔬 **Species Name:** **{species_name}**\n" + \
                           f"🧬 **OTU ID:** `{otu_id}`\n\n" + \
                           f"*Found in {filename}*"
                
                # If no exact match, try partial match (first 50 characters)
                partial_otu = otu_id[:50]
                partial_matches = df[df['otu_id'].str.upper().str.startswith(partial_otu)]
                
                if len(partial_matches) > 0:
                    species_name = partial_matches.iloc[0]['species_name']
                    return f"**OTU ID Lookup Result (Partial Match):**\n\n" + \
                           f"🧬 **Query OTU ID:** `{otu_id[:50]}...`\n" + \
                           f"🔬 **Species Name:** **{species_name}**\n\n" + \
                           f"*Found partial match in {filename}*"
        
        return f"**OTU ID Not Found:**\n\n" + \
               f"🧬 **Query OTU ID:** `{otu_id[:50]}...`\n" + \
               f"❌ No matching species found in uploaded files.\n\n" + \
               f"Please ensure you have uploaded a species mapping CSV file with 'otu_id' and 'species_name' columns."
    
    def lookup_unique_name_for_species(self, species_name: str) -> str:
        """Look up OTU ID for a given species name"""
        # Clean up the species name - remove common prefixes
        species_name = species_name.strip().lower()
        species_name = re.sub(r'^(the\s+)?(species\s+)?(name\s+)?', '', species_name)
        species_name = species_name.strip()
        
        # Search through uploaded CSV files for species mapping
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            # Check if this looks like a species mapping file with unique_name column
            if 'unique_name' in df.columns and 'species_name' in df.columns:
                # Look for exact match (case insensitive)
                matching_rows = df[df['species_name'].str.lower() == species_name]
                
                if len(matching_rows) > 0:
                    unique_name = matching_rows.iloc[0]['unique_name']
                    actual_species = matching_rows.iloc[0]['species_name']
                    
                    return f"**Unique Name Lookup Result:**\n\n" + \
                           f"🔬 **Species Name:** **{actual_species}**\n" + \
                           f"🏷️ **Unique Name:** **{unique_name}**\n\n" + \
                           f"*Found in {filename}*"
                
                # Try partial match
                partial_matches = df[df['species_name'].str.lower().str.contains(species_name, na=False)]
                
                if len(partial_matches) > 0:
                    results = []
                    for _, row in partial_matches.head(3).iterrows():  # Show top 3 matches
                        results.append(f"• **{row['species_name']}**: {row['unique_name']}")
                    
                    return f"**Partial Species Matches Found:**\n\n" + \
                           f"🔍 **Search term:** {species_name}\n\n" + \
                           "\n".join(results) + \
                           f"\n\n*Found in {filename}*"
        
        return f"**Species Not Found:**\n\n" + \
               f"🔍 **Search term:** {species_name}\n" + \
               f"❌ No matching unique name found in uploaded files.\n\n" + \
               f"Please ensure you have uploaded a species mapping CSV file with 'unique_name' and 'species_name' columns."
    
    def lookup_species_to_otu(self, species_name: str) -> str:
        """Look up OTU ID for a given species name"""
        # Clean up the species name - remove common prefixes
        species_name = species_name.strip().lower()
        species_name = re.sub(r'^(the\s+)?(species\s+)?(name\s+)?', '', species_name)
        species_name = species_name.strip()
        
        # Search through uploaded CSV files for species mapping
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            # Check if this looks like a species mapping file
            if 'otu_id' in df.columns and 'species_name' in df.columns:
                # Look for exact match (case insensitive)
                matching_rows = df[df['species_name'].str.lower() == species_name]
                
                if len(matching_rows) > 0:
                    otu_id = matching_rows.iloc[0]['otu_id']
                    actual_species = matching_rows.iloc[0]['species_name']
                    
                    # Include additional info if available
                    additional_info = ""
                    if 'unique_name' in df.columns:
                        unique_name = matching_rows.iloc[0]['unique_name']
                        additional_info += f"🏷️ **Unique Name:** {unique_name}\n"
                    
                    return f"**OTU ID Lookup Result:**\n\n" + \
                           f"🔬 **Species Name:** **{species_name}**\n" + \
                           f"🧬 **OTU ID:** `{otu_id}`\n" + \
                           additional_info + \
                           f"\n*Found in {filename}*"
                
                # Try partial match
                partial_matches = df[df['species_name'].str.lower().str.contains(species_name, na=False)]
                
                if len(partial_matches) > 0:
                    results = []
                    for _, row in partial_matches.head(3).iterrows():  # Show top 3 matches
                        results.append(f"• **{row['species_name']}**: `{row['otu_id']}`")
                    
                    return f"**Partial Species Matches Found:**\n\n" + \
                           f"🔍 **Search term:** {species_name}\n\n" + \
                           "\n".join(results) + \
                           f"\n\n*Found in {filename}*"
        
        return f"**Species Not Found:**\n\n" + \
               f"🔍 **Search term:** {species_name}\n" + \
               f"❌ No matching OTU ID found in uploaded files.\n\n" + \
               f"Please ensure you have uploaded a species mapping CSV file with 'otu_id' and 'species_name' columns."
    
    def lookup_numeric_id_for_species(self, species_name: str) -> str:
        """Look up numeric ID for a given species name"""
        import re
        
        # Clean up the species name - remove common prefixes
        species_name = species_name.strip().lower()
        species_name = re.sub(r'^(the\s+)?(species\s+)?(name\s+)?', '', species_name)
        species_name = species_name.strip()
        
        # Search through uploaded CSV files for species mapping
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            # Check if this looks like a species mapping file with numeric_id column
            if 'numeric_id' in df.columns and 'species_name' in df.columns:
                # Look for exact match (case insensitive)
                matching_rows = df[df['species_name'].str.lower() == species_name]
                
                if len(matching_rows) > 0:
                    numeric_id = matching_rows.iloc[0]['numeric_id']
                    actual_species = matching_rows.iloc[0]['species_name']
                    
                    # Include additional info if available
                    additional_info = ""
                    if 'otu_id' in df.columns:
                        otu_id = matching_rows.iloc[0]['otu_id']
                        additional_info += f"🧬 **OTU ID:** `{otu_id[:50]}{'...' if len(otu_id) > 50 else ''}`\n"
                    if 'unique_name' in df.columns:
                        unique_name = matching_rows.iloc[0]['unique_name']
                        additional_info += f"🏷️ **Unique Name:** {unique_name}\n"
                    
                    return f"**Numeric ID Lookup Result:**\n\n" + \
                           f"🔬 **Species Name:** **{species_name}**\n" + \
                           f"🔢 **Numeric ID:** **{numeric_id}**\n" + \
                           additional_info + \
                           f"\n*Found in {filename}*"
                
                # Try partial match
                partial_matches = df[df['species_name'].str.lower().str.contains(species_name, na=False)]
                
                if len(partial_matches) > 0:
                    results = []
                    for _, row in partial_matches.head(5).iterrows():  # Show top 5 matches
                        results.append(f"• **{row['species_name']}**: ID = {row['numeric_id']}")
                    
                    return f"**Partial Species Matches Found:**\n\n" + \
                           f"🔍 **Search term:** {species_name}\n\n" + \
                           "\n".join(results) + \
                           f"\n\n*Found in {filename}*"
        
        return f"**Species Not Found:**\n\n" + \
               f"🔍 **Search term:** {species_name}\n" + \
               f"❌ No matching numeric ID found in uploaded files.\n\n" + \
               f"Please ensure you have uploaded a species mapping CSV file with 'numeric_id' and 'species_name' columns."
    
    def lookup_taxonomy_info(self, search_term: str, query_lower: str) -> str:
        """Look up taxonomy information for a given genus/family/etc."""
        import re
        
        # Clean up the search term
        search_term = search_term.strip().lower()
        
        # Determine what type of lookup this is based on the query structure
        if 'family name of genus' in query_lower:
            return self.lookup_family_of_genus(search_term)
        elif 'genus of' in query_lower and 'family' in query_lower:
            return self.lookup_genus_of_family(search_term)
        elif 'phylum of order' in query_lower:
            return self.lookup_phylum_of_order(search_term)
        elif 'phylum of genus' in query_lower:
            return self.lookup_phylum_of_genus(search_term)
        elif 'class of genus' in query_lower:
            return self.lookup_class_of_genus(search_term)
        elif 'order of genus' in query_lower:
            return self.lookup_order_of_genus(search_term)
        elif 'order of phylum' in query_lower:
            return self.lookup_order_of_phylum(search_term)
        elif 'kingdom of' in query_lower:
            return self.lookup_kingdom_of_taxon(search_term)
        elif 'species of genus' in query_lower:
            return self.lookup_species_of_genus(search_term)
        elif 'species' in query_lower and any(rank in query_lower for rank in ['order', 'family', 'class', 'phylum', 'kingdom']):
            return self.lookup_species_by_taxonomy(search_term, query_lower)
        else:
            return self.lookup_family_of_genus(search_term)
    
    def lookup_family_of_genus(self, genus_name: str) -> str:
        """Look up family name for a given genus"""
        genus_name = genus_name.strip().lower()
        
        # Search through uploaded CSV files for taxonomy data
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            # Check if this looks like a taxonomy file
            if 'Genus' in df.columns and 'Family' in df.columns:
                # Look for exact match (case insensitive)
                matching_rows = df[df['Genus'].str.lower() == genus_name]
                
                if len(matching_rows) > 0:
                    family_name = matching_rows.iloc[0]['Family']
                    actual_genus = matching_rows.iloc[0]['Genus']
                    
                    # Include additional taxonomy info if available
                    additional_info = ""
                    if 'Phylum' in df.columns:
                        phylum = matching_rows.iloc[0]['Phylum']
                        additional_info += f"🌿 **Phylum:** {phylum}\n"
                    if 'Class' in df.columns:
                        class_name = matching_rows.iloc[0]['Class']
                        additional_info += f"📚 **Class:** {class_name}\n"
                    if 'Order' in df.columns:
                        order_name = matching_rows.iloc[0]['Order']
                        additional_info += f"📋 **Order:** {order_name}\n"
                    
                    return f"**Taxonomy Lookup Result:**\n\n" + \
                           f"🔬 **Genus:** **{actual_genus}**\n" + \
                           f"👨‍👩‍👧‍👦 **Family:** **{family_name}**\n" + \
                           additional_info + \
                           f"\n*Found in {filename}*"
                
                # Try partial match
                partial_matches = df[df['Genus'].str.lower().str.contains(genus_name, na=False)]
                
                if len(partial_matches) > 0:
                    results = []
                    for _, row in partial_matches.head(3).iterrows():
                        results.append(f"• **{row['Genus']}**: {row['Family']}")
                    
                    return f"**Partial Genus Matches Found:**\n\n" + \
                           f"🔍 **Search term:** {genus_name}\n\n" + \
                           "\n".join(results) + \
                           f"\n\n*Found in {filename}*"
        
        return f"**Genus Not Found:**\n\n" + \
               f"🔍 **Search term:** {genus_name}\n" + \
               f"❌ No matching genus found in uploaded taxonomy files.\n\n" + \
               f"Please ensure you have uploaded a taxonomy CSV file with 'Genus' and 'Family' columns."
    
    def lookup_genus_of_family(self, family_name: str) -> str:
        """Look up genus names for a given family"""
        family_name = family_name.strip().lower()
        
        # Search through uploaded CSV files for taxonomy data
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            # Check if this looks like a taxonomy file
            if 'Genus' in df.columns and 'Family' in df.columns:
                # Look for exact match (case insensitive)
                matching_rows = df[df['Family'].str.lower() == family_name]
                
                if len(matching_rows) > 0:
                    genera = matching_rows['Genus'].unique()
                    actual_family = matching_rows.iloc[0]['Family']
                    
                    result = f"**Taxonomy Lookup Result:**\n\n" + \
                             f"👨‍👩‍👧‍👦 **Family:** **{actual_family}**\n" + \
                             f"🔬 **Genera ({len(genera)}):**\n\n"
                    
                    for genus in sorted(genera):
                        result += f"• {genus}\n"
                    
                    result += f"\n*Found in {filename}*"
                    return result
                
                # Try partial match
                partial_matches = df[df['Family'].str.lower().str.contains(family_name, na=False)]
                
                if len(partial_matches) > 0:
                    families = partial_matches['Family'].unique()[:3]
                    results = []
                    for family in families:
                        genera_in_family = partial_matches[partial_matches['Family'] == family]['Genus'].unique()
                        results.append(f"• **{family}**: {', '.join(genera_in_family[:3])}{'...' if len(genera_in_family) > 3 else ''}")
                    
                    return f"**Partial Family Matches Found:**\n\n" + \
                           f"🔍 **Search term:** {family_name}\n\n" + \
                           "\n".join(results) + \
                           f"\n\n*Found in {filename}*"
        
        return f"**Family Not Found:**\n\n" + \
               f"🔍 **Search term:** {family_name}\n" + \
               f"❌ No matching family found in uploaded taxonomy files.\n\n" + \
               f"Please ensure you have uploaded a taxonomy CSV file with 'Genus' and 'Family' columns."
    
    def lookup_phylum_of_genus(self, genus_name: str) -> str:
        """Look up phylum for a given genus"""
        genus_name = genus_name.strip().lower()
        
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            if 'Genus' in df.columns and 'Phylum' in df.columns:
                matching_rows = df[df['Genus'].str.lower() == genus_name]
                
                if len(matching_rows) > 0:
                    phylum_name = matching_rows.iloc[0]['Phylum']
                    actual_genus = matching_rows.iloc[0]['Genus']
                    
                    return f"**Taxonomy Lookup Result:**\n\n" + \
                           f"🔬 **Genus:** **{actual_genus}**\n" + \
                           f"🌿 **Phylum:** **{phylum_name}**\n" + \
                           f"\n*Found in {filename}*"
        
        return f"**Genus Not Found:**\n\n" + \
               f"🔍 **Search term:** {genus_name}\n" + \
               f"❌ No matching genus found in uploaded taxonomy files."
    
    def lookup_class_of_genus(self, genus_name: str) -> str:
        """Look up class for a given genus"""
        genus_name = genus_name.strip().lower()
        
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            if 'Genus' in df.columns and 'Class' in df.columns:
                matching_rows = df[df['Genus'].str.lower() == genus_name]
                
                if len(matching_rows) > 0:
                    class_name = matching_rows.iloc[0]['Class']
                    actual_genus = matching_rows.iloc[0]['Genus']
                    
                    return f"**Taxonomy Lookup Result:**\n\n" + \
                           f"🔬 **Genus:** **{actual_genus}**\n" + \
                           f"📚 **Class:** **{class_name}**\n" + \
                           f"\n*Found in {filename}*"
        
        return f"**Genus Not Found:**\n\n" + \
               f"🔍 **Search term:** {genus_name}\n" + \
               f"❌ No matching genus found in uploaded taxonomy files."
    
    def lookup_order_of_genus(self, genus_name: str) -> str:
        """Look up order for a given genus"""
        genus_name = genus_name.strip().lower()
        
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            if 'Genus' in df.columns and 'Order' in df.columns:
                matching_rows = df[df['Genus'].str.lower() == genus_name]
                
                if len(matching_rows) > 0:
                    order_name = matching_rows.iloc[0]['Order']
                    actual_genus = matching_rows.iloc[0]['Genus']
                    
                    return f"**Taxonomy Lookup Result:**\n\n" + \
                           f"🔬 **Genus:** **{actual_genus}**\n" + \
                           f"📋 **Order:** **{order_name}**\n" + \
                           f"\n*Found in {filename}*"
        
        return f"**Genus Not Found:**\n\n" + \
               f"🔍 **Search term:** {genus_name}\n" + \
               f"❌ No matching genus found in uploaded taxonomy files."
    
    def lookup_phylum_of_order(self, order_name: str) -> str:
        """Look up phylum for a given order"""
        order_name = order_name.strip().lower()
        
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            if 'Order' in df.columns and 'Phylum' in df.columns:
                matching_rows = df[df['Order'].str.lower() == order_name]
                
                if len(matching_rows) > 0:
                    phylum_name = matching_rows.iloc[0]['Phylum']
                    actual_order = matching_rows.iloc[0]['Order']
                    
                    # Include additional taxonomy info if available
                    additional_info = ""
                    if 'Class' in df.columns:
                        class_name = matching_rows.iloc[0]['Class']
                        additional_info += f"📚 **Class:** {class_name}\n"
                    
                    return f"**Taxonomy Lookup Result:**\n\n" + \
                           f"📋 **Order:** **{actual_order}**\n" + \
                           f"🌿 **Phylum:** **{phylum_name}**\n" + \
                           additional_info + \
                           f"\n*Found in {filename}*"
                
                # Try partial match
                partial_matches = df[df['Order'].str.lower().str.contains(order_name, na=False)]
                
                if len(partial_matches) > 0:
                    results = []
                    for _, row in partial_matches.head(3).iterrows():
                        results.append(f"• **{row['Order']}**: {row['Phylum']}")
                    
                    return f"**Partial Order Matches Found:**\n\n" + \
                           f"🔍 **Search term:** {order_name}\n\n" + \
                           "\n".join(results) + \
                           f"\n\n*Found in {filename}*"
        
        return f"**Order Not Found:**\n\n" + \
               f"🔍 **Search term:** {order_name}\n" + \
               f"❌ No matching order found in uploaded taxonomy files.\n\n" + \
               f"Please ensure you have uploaded a taxonomy CSV file with 'Order' and 'Phylum' columns."
    
    def lookup_kingdom_of_taxon(self, taxon_name: str) -> str:
        """Look up kingdom for a given taxon (genus, species, etc.)"""
        taxon_name = taxon_name.strip().lower()
        
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            if 'Kingdom' in df.columns:
                # Try matching against different taxonomic levels
                for col in ['Genus', 'Species', 'Family', 'Order', 'Class']:
                    if col in df.columns:
                        matching_rows = df[df[col].str.lower() == taxon_name]
                        
                        if len(matching_rows) > 0:
                            kingdom_name = matching_rows.iloc[0]['Kingdom']
                            actual_taxon = matching_rows.iloc[0][col]
                            
                            return f"**Taxonomy Lookup Result:**\n\n" + \
                                   f"🏛️ **{col}:** **{actual_taxon}**\n" + \
                                   f"👑 **Kingdom:** **{kingdom_name}**\n" + \
                                   f"\n*Found in {filename}*"
        
        return f"**Taxon Not Found:**\n\n" + \
               f"🔍 **Search term:** {taxon_name}\n" + \
               f"❌ No matching taxon found in uploaded taxonomy files."
    
    def lookup_species_of_genus(self, genus_name: str) -> str:
        """Look up species for a given genus"""
        genus_name = genus_name.strip().lower()
        
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            if 'Genus' in df.columns and 'Species' in df.columns:
                matching_rows = df[df['Genus'].str.lower() == genus_name]
                
                if len(matching_rows) > 0:
                    species_list = matching_rows['Species'].unique()
                    actual_genus = matching_rows.iloc[0]['Genus']
                    
                    result = f"**Taxonomy Lookup Result:**\n\n" + \
                             f"🔬 **Genus:** **{actual_genus}**\n" + \
                             f"🧬 **Species ({len(species_list)}):**\n\n"
                    
                    for species in sorted(species_list):
                        if pd.notna(species) and species.strip():
                            result += f"• {species}\n"
                    
                    result += f"\n*Found in {filename}*"
                    return result
        
        return f"**Genus Not Found:**\n\n" + \
               f"🔍 **Search term:** {genus_name}\n" + \
               f"❌ No matching genus found in uploaded taxonomy files."
    
    def lookup_species_by_taxonomy(self, taxon_name: str, query_lower: str) -> str:
        """Look up species that belong to a specific taxonomic group"""
        taxon_name = taxon_name.strip().lower()
        
        # Determine which taxonomic rank we're searching by
        rank_column = None
        if 'order' in query_lower:
            rank_column = 'Order'
        elif 'family' in query_lower:
            rank_column = 'Family'
        elif 'class' in query_lower:
            rank_column = 'Class'
        elif 'phylum' in query_lower:
            rank_column = 'Phylum'
        elif 'kingdom' in query_lower:
            rank_column = 'Kingdom'
        
        if not rank_column:
            return "Could not determine taxonomic rank from query."
        
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            if rank_column in df.columns and 'Species' in df.columns:
                matching_rows = df[df[rank_column].str.lower() == taxon_name]
                
                if len(matching_rows) > 0:
                    species_list = matching_rows['Species'].unique()
                    actual_taxon = matching_rows.iloc[0][rank_column]
                    
                    # Filter out empty/null species
                    valid_species = [s for s in species_list if pd.notna(s) and s.strip()]
                    
                    result = f"**Species Lookup Result:**\n\n" + \
                             f"📋 **{rank_column}:** **{actual_taxon}**\n" + \
                             f"🧬 **Species ({len(valid_species)}):**\n\n"
                    
                    for species in sorted(valid_species):
                        result += f"• {species}\n"
                    
                    result += f"\n*Found in {filename}*"
                    return result
                
                # Try partial match
                partial_matches = df[df[rank_column].str.lower().str.contains(taxon_name, na=False)]
                
                if len(partial_matches) > 0:
                    taxa = partial_matches[rank_column].unique()[:3]
                    results = []
                    for taxon in taxa:
                        species_count = len(partial_matches[partial_matches[rank_column] == taxon]['Species'].unique())
                        results.append(f"• **{taxon}**: {species_count} species")
                    
                    return f"**Partial {rank_column} Matches Found:**\n\n" + \
                           f"🔍 **Search term:** {taxon_name}\n\n" + \
                           "\n".join(results) + \
                           f"\n\n*Found in {filename}*"
        
        return f"**{rank_column} Not Found:**\n\n" + \
               f"🔍 **Search term:** {taxon_name}\n" + \
               f"❌ No matching {rank_column.lower()} found in uploaded taxonomy files.\n\n" + \
               f"Please ensure you have uploaded a taxonomy CSV file with '{rank_column}' and 'Species' columns."
    
    def lookup_taxonomy_by_otu_sequence(self, otu_sequence: str, rank_requested: str) -> str:
        """Look up taxonomy information for an OTU sequence"""
        # First, find the species name for this OTU sequence
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            # Check species mapping files first
            if 'otu_id' in df.columns and 'species_name' in df.columns:
                # Try exact match
                matching_rows = df[df['otu_id'].str.upper() == otu_sequence]
                
                if len(matching_rows) > 0:
                    species_name = matching_rows.iloc[0]['species_name']
                    
                    # Now look up taxonomy info for this species
                    return self.lookup_taxonomy_for_species(species_name, rank_requested, otu_sequence)
                
                # Try partial match (first 100 characters)
                partial_otu = otu_sequence[:100]
                partial_matches = df[df['otu_id'].str.upper().str.startswith(partial_otu)]
                
                if len(partial_matches) > 0:
                    species_name = partial_matches.iloc[0]['species_name']
                    return self.lookup_taxonomy_for_species(species_name, rank_requested, otu_sequence)
        
        return f"**OTU Sequence Not Found:**\n\n" + \
               f"🧬 **Query OTU:** `{otu_sequence[:50]}...`\n" + \
               f"❌ No matching OTU sequence found in uploaded files.\n\n" + \
               f"Please ensure you have uploaded a species mapping CSV file with 'otu_id' and 'species_name' columns."
    
    def lookup_taxonomy_for_species(self, species_name: str, rank_requested: str, otu_sequence: str) -> str:
        """Look up specific taxonomy rank for a species"""
        rank_column_map = {
            'order': 'Order',
            'family': 'Family', 
            'class': 'Class',
            'phylum': 'Phylum',
            'kingdom': 'Kingdom'
        }
        
        rank_column = rank_column_map.get(rank_requested.lower())
        if not rank_column:
            return f"Unknown taxonomy rank: {rank_requested}"
        
        # Search taxonomy files for this species
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            # Check if this is a taxonomy file with the species and requested rank
            if 'Species' in df.columns and rank_column in df.columns:
                # Try to match by species name
                matching_rows = df[df['Species'].str.lower().str.contains(species_name.lower(), na=False)]
                
                if len(matching_rows) > 0:
                    rank_value = matching_rows.iloc[0][rank_column]
                    actual_species = matching_rows.iloc[0]['Species']
                    
                    return f"**OTU Taxonomy Lookup Result:**\n\n" + \
                           f"🧬 **OTU Sequence:** `{otu_sequence[:50]}...`\n" + \
                           f"🔬 **Species:** **{species_name}**\n" + \
                           f"📋 **{rank_column}:** **{rank_value}**\n" + \
                           f"\n*Found in {filename}*"
        
        return f"**Taxonomy Not Found:**\n\n" + \
               f"🔬 **Species:** {species_name}\n" + \
               f"❌ No taxonomy information found for this species.\n\n" + \
               f"Please ensure you have uploaded a taxonomy CSV file with 'Species' and '{rank_column}' columns."
    
    def lookup_order_of_phylum(self, phylum_name: str) -> str:
        """Look up orders for a given phylum"""
        phylum_name = phylum_name.strip().lower()
        
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']
            
            if 'Phylum' in df.columns and 'Order' in df.columns:
                matching_rows = df[df['Phylum'].str.lower() == phylum_name]
                
                if len(matching_rows) > 0:
                    orders = matching_rows['Order'].unique()
                    actual_phylum = matching_rows.iloc[0]['Phylum']
                    
                    result = f"**Taxonomy Lookup Result:**\n\n" + \
                             f"🌿 **Phylum:** **{actual_phylum}**\n" + \
                             f"📋 **Orders ({len(orders)}):**\n\n"
                    
                    for order in sorted(orders):
                        if pd.notna(order) and order.strip():
                            result += f"• {order}\n"
                    
                    result += f"\n*Found in {filename}*"
                    return result
        
        return f"**Phylum Not Found:**\n\n" + \
               f"🔍 **Search term:** {phylum_name}\n" + \
               f"❌ No matching phylum found in uploaded taxonomy files."
    
    def find_most_relevant_csv(self, query_lower: str):
        """Find the most relevant CSV file for a given query"""
        if not self.uploaded_files['csv_files']:
            return None
        
        # Return the most recently uploaded CSV file based on timestamp
        latest_file = max(self.uploaded_files['csv_files'].items(), 
                         key=lambda x: x[1]['timestamp'])
        return latest_file
    
    def handle_intelligent_data_query(self, df, filename: str, csv_type: str, query_lower: str) -> str:
        """Intelligently handle any query about CSV data"""
        import re
        
        # Extract column names and data types
        columns = df.columns.tolist()
        numeric_columns = df.select_dtypes(include=['number']).columns.tolist()
        text_columns = df.select_dtypes(include=['object']).columns.tolist()
        
        # 1. Handle specific column queries (but avoid generic column names in species queries)
        for col in columns:
            if col.lower() in query_lower and not any(term in query_lower for term in ['unique_name', 'otu_id', 'species']):
                return self.analyze_specific_column(df, col, filename, query_lower)
            
        # 2. Handle "unique" queries for any column
        if 'unique' in query_lower:
            # Try to find which column the user is asking about
            for col in text_columns:
                if any(word in query_lower for word in col.lower().split('_')):
                    return self.get_unique_values(df, col, filename)
                
                # If no specific column found, show unique values for all text columns
                return self.get_all_unique_values(df, filename)
            
        # 3. Handle "count" or "how many" queries
        if any(term in query_lower for term in ['count', 'how many', 'number of']):
            return self.handle_count_queries(df, filename, query_lower)
            
        # 4. Handle "show", "list", "display" queries
        if any(term in query_lower for term in ['show', 'list', 'display', 'what']):
            return self.handle_display_queries(df, filename, query_lower)
            
        # 5. Handle statistical queries
        if any(term in query_lower for term in ['average', 'mean', 'max', 'min', 'sum', 'statistics']):
            return self.handle_statistical_queries(df, filename, query_lower, numeric_columns)
            
        # 6. Handle search queries (find, search, contains)
        if any(term in query_lower for term in ['find', 'search', 'contains', 'where']):
            # Extract search term
            search_terms = query_lower.replace('find', '').replace('search', '').strip()
            if search_terms:
                return self.handle_search_queries(df, filename, search_terms)
            else:
                return "Please specify what you want to search for. Example: 'find Bacillus'"
            
        # 7. Handle correlation and relationship queries
        if any(term in query_lower for term in ['correlation', 'relationship', 'related']):
            return self.handle_correlation_queries(df, filename, numeric_columns)
            
        # 8. Default: provide data overview
        return self.get_data_overview(df, filename)
            
    def analyze_specific_column(self, df, column: str, filename: str, query_lower: str) -> str:
        """Analyze a specific column based on the query"""
        col_data = df[column].dropna()
        
        if col_data.empty:
            return f"**Column '{column}' in {filename}:**\n\n❌ No data available (all values are null)"
        
        result = f"**Analysis of '{column}' in {filename}:**\n\n"
        
        if col_data.dtype in ['object', 'string']:
            # Text column analysis
            unique_count = col_data.nunique()
            total_count = len(col_data)
            
            result += f"📊 **Total values:** {total_count}\n"
            result += f"🔢 **Unique values:** {unique_count}\n\n"
            
            if unique_count <= 20:
                result += "**All unique values:**\n"
                for value in sorted(col_data.unique()):
                    count = (col_data == value).sum()
                    result += f"• {value} ({count} occurrences)\n"
            else:
                result += "**Top 10 most common values:**\n"
                value_counts = col_data.value_counts().head(10)
                for value, count in value_counts.items():
                    result += f"• {value} ({count} occurrences)\n"
        else:
            # Numeric column analysis
            result += f"📊 **Count:** {len(col_data)}\n"
            result += f"📈 **Mean:** {col_data.mean():.2f}\n"
            result += f"📉 **Min:** {col_data.min()}\n"
            result += f"📊 **Max:** {col_data.max()}\n"
            result += f"🎯 **Median:** {col_data.median():.2f}\n"
            result += f"📏 **Std Dev:** {col_data.std():.2f}\n"
        
        return result
    
    def get_unique_values(self, df, column: str, filename: str) -> str:
        """Get unique values for a specific column"""
        if column not in df.columns:
            return f"Column '{column}' not found in {filename}. Available columns: {', '.join(df.columns)}"
        
        unique_values = df[column].dropna().unique()
        
        result = f"**Unique values in '{column}' ({filename}):**\n\n"
        result += f"📊 **Total unique values:** {len(unique_values)}\n\n"
        
        if len(unique_values) <= 50:
            result += "**All values:**\n"
            for value in sorted(unique_values):
                result += f"• {value}\n"
        else:
            result += "**First 50 values:**\n"
            for value in sorted(unique_values)[:50]:
                result += f"• {value}\n"
            result += f"\n*... and {len(unique_values) - 50} more values*"
        
        return result
    
    def get_all_unique_values(self, df, filename: str) -> str:
        """Get unique values for all text columns"""
        text_columns = df.select_dtypes(include=['object']).columns.tolist()
        
        if not text_columns:
            return f"No text columns found in {filename} for unique value analysis."
        
        result = f"**Unique Values Summary for {filename}:**\n\n"
        
        for col in text_columns[:5]:  # Limit to first 5 columns
            unique_count = df[col].nunique()
            result += f"📊 **{col}:** {unique_count} unique values\n"
        
        return result
    
    def handle_count_queries(self, df, filename: str, query_lower: str) -> str:
        """Handle counting queries"""
        result = f"**Count Information for {filename}:**\n\n"
        result += f"📋 **Total rows:** {len(df)}\n"
        result += f"📊 **Total columns:** {len(df.columns)}\n\n"
        
        # Count non-null values per column
        result += "**Non-null counts per column:**\n"
        for col in df.columns[:10]:  # Show first 10 columns
            non_null_count = df[col].count()
            result += f"• {col}: {non_null_count}\n"
        
        if len(df.columns) > 10:
            result += f"\n*... and {len(df.columns) - 10} more columns*"
        
        return result
    
    def handle_display_queries(self, df, filename: str, query_lower: str) -> str:
        """Handle display/show queries"""
        result = f"**Data Preview for {filename}:**\n\n"
        result += f"📊 **Shape:** {df.shape[0]} rows × {df.shape[1]} columns\n\n"
        
        result += "**Column Information:**\n"
        for col in df.columns:
            dtype = str(df[col].dtype)
            non_null = df[col].count()
            result += f"• **{col}** ({dtype}): {non_null} non-null values\n"
        
        result += "\n**Sample Data (first 3 rows):**\n"
        for i, (_, row) in enumerate(df.head(10).iterrows()):
            result += f"\n**Row {i+1}:**\n"
            for col in df.columns[:5]:  # Show first 5 columns
                value = str(row[col])[:50]  # Truncate long values
                result += f"  • {col}: {value}\n"
        
        return result
    
    def handle_statistical_queries(self, df, filename: str, query_lower: str, numeric_columns: list) -> str:
        """Handle statistical analysis queries"""
        if not numeric_columns:
            return f"No numeric columns found in {filename} for statistical analysis."
        
        result = f"**Statistical Summary for {filename}:**\n\n"
        
        for col in numeric_columns[:5]:  # Limit to first 5 numeric columns
            col_data = df[col].dropna()
            if len(col_data) > 0:
                result += f"**{col}:**\n"
                result += f"  • Count: {len(col_data)}\n"
                result += f"  • Mean: {col_data.mean():.2f}\n"
                result += f"  • Std: {col_data.std():.2f}\n"
                result += f"  • Min: {col_data.min()}\n"
                result += f"  • Max: {col_data.max()}\n\n"
        
        return result
    
    def handle_search_queries(self, df, filename: str, query_lower: str) -> str:
        """Handle search and filter queries"""
        # Extract search terms from the query
        import re
        
        # Look for quoted search terms
        quoted_terms = re.findall(r'["\']([^"\']*)["\']', query_lower)
        
        if not quoted_terms:
            return "**Search Help for " + filename + ":**\\n\\nTo search the data, please specify search terms in quotes.\\nExample: 'find \"microvirga\"' or 'search \"bacteria\"'"
        
        search_term = quoted_terms[0]
        results = []
        
        # Search in all text columns
        text_columns = df.select_dtypes(include=['object']).columns
        for col in text_columns:
            matches = df[df[col].astype(str).str.contains(search_term, case=False, na=False)]
            if not matches.empty:
                results.append(f"**{col}:** {len(matches)} matches")
                # Show first few matches
                for idx, row in matches.head(5).iterrows():
                    results.append(f"  • {row[col]}")
                if len(matches) > 5:
                    results.append(f"  • *...and {len(matches) - 5} more matches*")
        
        if results:
            return f"**Search Results for '{search_term}' in {filename}:**\n\n" + "\n".join(results)
        else:
            return f"No matches found for '{search_term}' in {filename}."
    
    def handle_correlation_queries(self, df, filename: str, numeric_columns: list) -> str:
        """Handle correlation analysis queries"""
        if len(numeric_columns) < 2:
            return f"Need at least 2 numeric columns for correlation analysis. Found {len(numeric_columns)} in {filename}."
        
        # Calculate correlations
        corr_matrix = df[numeric_columns].corr()
        
        result = f"**Correlation Analysis for {filename}:**\n\n"
        
        # Find strongest correlations
        strong_corrs = []
        for i in range(len(numeric_columns)):
            for j in range(i+1, len(numeric_columns)):
                corr_val = corr_matrix.iloc[i, j]
                if abs(corr_val) > 0.5:  # Strong correlation threshold
                    strong_corrs.append((numeric_columns[i], numeric_columns[j], corr_val))
        
        if strong_corrs:
            result += "**Strong correlations (|r| > 0.5):**\n"
            for col1, col2, corr in sorted(strong_corrs, key=lambda x: abs(x[2]), reverse=True):
                result += f"• {col1} ↔ {col2}: {corr:.3f}\n"
        else:
            result += "No strong correlations found (all |r| ≤ 0.5)\n"
        
        return result
    
    def get_data_overview(self, df, filename: str) -> str:
        """Provide a comprehensive data overview"""
        result = f"**Data Overview for {filename}:**\n\n"
        result += f"📊 **Shape:** {df.shape[0]} rows × {df.shape[1]} columns\n\n"
        
        # Column types
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        text_cols = df.select_dtypes(include=['object']).columns.tolist()
        
        result += f"**Column Types:**\n"
        result += f"• Numeric columns: {len(numeric_cols)}\n"
        result += f"• Text columns: {len(text_cols)}\n\n"
        
        result += "**Available Data:**\n"
        result += "You can ask questions like:\n"
        result += "• 'Show unique values in [column_name]'\n"
        result += "• 'What is the average [numeric_column]?'\n"
        result += "• 'Find rows containing \"search_term\"'\n"
        result += "• 'How many records are there?'\n"
        result += "• 'Show correlation between columns'\n"
        
        return result

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
    response = chatbot.process_user_query(user_message)
    
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
        filename = data.get('filename', 'uploaded_data.csv')
        
        # Store the CSV data in uploaded_files for query handling
        chatbot.uploaded_files['csv_files'][filename] = {
            'path': filename,
            'data': otu_data,
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        
        result = chatbot.analyze_microbiome_data(otu_data, depth)
        result['file_stored'] = True
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

@app.route('/upload_depth', methods=['POST'])
def upload_depth():
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': "No file part named 'file'"}), 400
        f = request.files['file']
        if not f.filename:
            return jsonify({'success': False, 'error': "No file selected"}), 400

        os.makedirs(os.path.dirname(DEPTH_CSV_DEFAULT), exist_ok=True)
        f.save(DEPTH_CSV_DEFAULT)

        chatbot.load_depth_csv(DEPTH_CSV_DEFAULT)

        if not hasattr(chatbot, 'uploaded_files'):
            chatbot.uploaded_files = {'csv_files': {}, 'pdf_files': {}}

        # ✅ register the uploaded file in-memory so depth queries work immediately
        chatbot.uploaded_files['csv_files']['depth.csv'] = {
            'path': DEPTH_CSV_DEFAULT,
            'data': chatbot.depth_df.copy(),
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }

        return jsonify({'success': True, 'schema': chatbot.depth_schema})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 400

@app.route('/upload_file', methods=['POST'])
def upload_file():
    """Handle single file upload (PDF or CSV)"""
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
        
        if file_ext == 'pdf':
            # Process PDF file
            result = chatbot.process_pdf_file(filename, file_path)
        elif file_ext in ['csv', 'tsv', 'txt']:
            # Process CSV file
            result = chatbot.analyze_csv_file(filename, file_path)
        else:
            return jsonify({'success': False, 'error': 'Unsupported file type'}), 400
        
        # Add to chat history if successful
        if result.get('success'):
            chatbot.chat_history.append({
                'user': f'Uploaded {filename}',
                'bot': result['message'],
                'timestamp': datetime.now().isoformat()
            })
        
        return jsonify(result)
        
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
            return jsonify({'success': False, 'error': "Query is required"}), 400
        
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
