import os
import json
import numpy as np
import difflib
import pandas as pd
from flask import Flask, render_template, request, jsonify
import subprocess
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from datetime import datetime
import re
from typing import Dict, List, Tuple, Any
import signal
import sys
import PyPDF2
import fitz  
from werkzeug.utils import secure_filename

# Enable automatic conversion between pandas and R

app = Flask(__name__)
app.secret_key = os.environ.get("FLASK_SECRET_KEY", "dev-secret")

DEPTH_CSV_DEFAULT = os.environ.get("DEPTH_CSV", os.path.join("data","depth.csv"))

class MicrobiomeXAIChatbot:
    def _suggest_closest_taxa(self, df, column, name, max_suggestions=5):
        """Return a 'Did you mean' string using difflib if close matches exist."""
        try:
            vals = [str(v) for v in df[column].dropna().unique() if str(v).strip()]
            matches = difflib.get_close_matches(name, vals, n=max_suggestions, cutoff=0.6)
            if matches:
                return " Did you mean: " + ", ".join(matches) + "?"
        except Exception:
            pass
        return ""

    def _load_csv_from_known_locations(self, name_hint):
        """Try to find a CSV by name in uploaded files or self.data_dir."""
        # To look in uploaded files first
        for fname, info in self.uploaded_files.get('csv_files', {}).items():
            if name_hint.lower() in fname.lower():
                return info.get('data'), fname
        # To look in data directory
        roots = [getattr(self, 'data_dir', 'data'), '.']
        for root in roots:
            try:
                for fn in os.listdir(root):
                    if fn.lower().endswith('.csv') and name_hint.lower() in fn.lower():
                        try:
                            return pd.read_csv(os.path.join(root, fn)), fn
                        except Exception:
                            continue
            except Exception:
                continue
        return None, None

    def compare_interactions(self, a_df, b_df):
        """Return overlap/changes between two interaction tables (columns sp1, sp2, lnk optionally)."""
        def norm(df):
            df = df.copy()
            if 'sp1' not in df.columns or 'sp2' not in df.columns:
                rename_map = {}
                if 'species1' in df.columns: rename_map['species1'] = 'sp1'
                if 'species2' in df.columns: rename_map['species2'] = 'sp2'
                if 'interaction' in df.columns and 'lnk' not in df.columns: rename_map['interaction'] = 'lnk'
                df = df.rename(columns=rename_map)
            if 'lnk' not in df.columns: df['lnk'] = 'unknown'
            undir = set(tuple(sorted((str(r['sp1']), str(r['sp2'])))) for _, r in df.iterrows())
            return df, undir

        A, A_undir = norm(a_df)
        B, B_undir = norm(b_df)

        overlap_edges = A_undir & B_undir
        added_edges = B_undir - A_undir
        removed_edges = A_undir - B_undir

        def labels_for(df):
            d = {}
            for _, r in df.iterrows():
                key = tuple(sorted((str(r.get('sp1','')), str(r.get('sp2','')))))
                d.setdefault(key, set()).add(str(r.get('lnk','unknown')))
            return d

        labA = labels_for(A)
        labB = labels_for(B)
        changed = {k: {'a': sorted(labA.get(k, [])), 'b': sorted(labB.get(k, []))}
                   for k in overlap_edges if labA.get(k, set()) != labB.get(k, set())}

        union = A_undir | B_undir
        jaccard = (len(overlap_edges) / len(union)) if union else 0.0

        return {
            'a_edges': len(A_undir),
            'b_edges': len(B_undir),
            'overlap_edges': len(overlap_edges),
            'jaccard': jaccard,
            'added_examples': list(sorted(list(added_edges)))[:20],
            'removed_examples': list(sorted(list(removed_edges)))[:20],
            'label_changes': changed
        }

    def compare_taxonomy(self, a_df, b_df):
        """Compare two taxonomy tables by parent->child pairs for adjacent levels."""
        cols = [c for c in a_df.columns if c in b_df.columns]
        cols = [c for c in cols if c.lower() in ['kingdom','phylum','class','order','family','genus','species']]
        if len(cols) < 2:
            return {'error': 'Taxonomy files need at least two common taxonomic columns to compare.'}
        def pairs(df, parent, child):
            df = df.dropna(subset=[parent, child]).copy()
            df[parent] = df[parent].astype(str).str.strip()
            df[child]  = df[child].astype(str).str.strip()
            return set(tuple(x) for x in df[[parent, child]].drop_duplicates().values.tolist())
        results = {'columns_compared': cols, 'levels': []}
        for i in range(len(cols)-1):
            parent, child = cols[i], cols[i+1]
            a_pairs = pairs(a_df, parent, child)
            b_pairs = pairs(b_df, parent, child)
            overlap = a_pairs & b_pairs
            union = a_pairs | b_pairs
            results['levels'].append({
                'relation': f'{parent}->{child}',
                'a_pairs': len(a_pairs),
                'b_pairs': len(b_pairs),
                'overlap': len(overlap),
                'jaccard': (len(overlap)/len(union) if union else 0.0),
                'added_examples': list(sorted(list(b_pairs - a_pairs)))[:20],
                'removed_examples': list(sorted(list(a_pairs - b_pairs)))[:20],
            })
        return results

    def compare_mappings(self, a_df, b_df):
        """Compare two species mapping tables by unique_name/species_name coverage and overlap."""
        def keyset(df):
            key = None
            for cand in ['unique_name','species_name']:
                if cand in df.columns:
                    key = cand
                    break
            if not key:
                return set(), None
            vals = set(str(v).strip() for v in df[key].dropna().unique())
            return vals, key
        a_keys, a_col = keyset(a_df)
        b_keys, b_col = keyset(b_df)
        if a_col is None or b_col is None:
            return {'error': 'Both mapping files need a unique_name or species_name column.'}
        union = a_keys | b_keys
        return {
            'a_keys': len(a_keys),
            'b_keys': len(b_keys),
            'overlap': len(a_keys & b_keys),
            'jaccard': (len(a_keys & b_keys)/len(union) if union else 0.0),
            'only_in_a': list(sorted(list(a_keys - b_keys)))[:20],
            'only_in_b': list(sorted(list(b_keys - a_keys)))[:20],
        }

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
        self.allowed_extensions = {'csv','tsv','txt','pdf'}
        self.upload_dir = 'uploads'
        os.makedirs(self.upload_dir, exist_ok=True)
        
    def find_r_executable(self):
        """Find R executable on Windows system"""
        possible_paths = [
            r"C:\Program Files\R\R-*\bin\R.exe",
            r"C:\Program Files\R\R-*\bin\x64\R.exe", 
            r"C:\Program Files (x86)\R\R-*\bin\R.exe",
            r"C:\Program Files (x86)\R\R-*\bin\x64\R.exe",
            "R.exe", 
            "R"       
        ]
        
        import glob
        for pattern in possible_paths:
            if '*' in pattern:
                matches = glob.glob(pattern)
                if matches:
                    
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
            
            # Tests if R is working
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
            # For Preprocess the data to handle missing values and ensure proper format
            processed_data = self.preprocess_otu_data(otu_data)
            if isinstance(processed_data, dict) and 'error' in processed_data:
                return processed_data
            
            # Creating R script for analysis with optimized parameters using unique temp files
            input_path = os.path.join(os.getcwd(), 'input.csv')
            processed_data.to_csv(input_path, index=True)
            
            interactions_output = os.path.join(os.getcwd(), 'interactions.csv')
            interactions_enhanced_output = os.path.join(os.getcwd(), 'interactions_enhanced.csv')
            
            # Fix Windows compatibility: mclapply doesn't support multiple cores on Windows
            if os.name == 'nt':  # Windows
                ncores = 1
            else:
                ncores = min(6, os.cpu_count() or 1)  # Unix/Linux systems
            
            # Get optimized nperms from preprocessing 
            nperms = getattr(processed_data, '_nperms', 10)
            
            depth_str = ', '.join(map(str, depth)) if depth else ''
            
            # Converted paths to forward slashes for R compatibility
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
            capture_output=True, text=True, timeout=timeout_seconds, cwd=os.getcwd()
            )
        
            # Better error analysis
            if result.returncode != 0:
                error_msg = result.stderr.strip()
            
                # Categorizes common errors
                if "InfIntE" in error_msg:
                    return {"error": "InfIntE package error. Please ensure InfIntE is properly installed in R."}
                elif "Python3" in error_msg or "python" in error_msg.lower():
                    return {"error": "Python dependency issue. Please install numpy, texttable, and cython for PyGol."}
                elif "memory" in error_msg.lower() or "cannot allocate" in error_msg.lower():
                    return {"error": f"Memory error. Dataset too large ({original_rows}x{original_cols}). Try reducing dataset size."}
                elif "timeout" in error_msg.lower():
                    return {"error": f"Analysis timed out after {timeout_seconds//60} minutes. Try with a smaller dataset."}
                else:
                    return {"error": f"R analysis failed: {error_msg[:200]}..."}
        
            
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
        except FileNotFoundError:
            return {"error": "R executable not found. Please install R and add it to PATH."}
        except PermissionError:
            return {"error": "Permission denied accessing R or temp executable. Please check file permissions."}
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
            # If first column looks like sample names/IDs and has unique values, uses it as index
            first_col = processed_data.columns[0]
            if (processed_data[first_col].dtype == 'object' and 
                processed_data[first_col].nunique() == len(processed_data) and
                len(processed_data) > 1):
                processed_data = processed_data.set_index(first_col)
            
            # Remove completely empty rows and columns
            processed_data = processed_data.dropna(how='all', axis=0)  # For Removing empty rows
            processed_data = processed_data.dropna(how='all', axis=1)  # For Removing empty columns
            
            # Check if we still have data after removing empty rows/columns
            if processed_data.empty:
                return {"error": "Dataset contains only empty values"}
            
            # Handle missing values by filling with 0 (common practice we should do for abundance data)
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
            if processed_data.shape[0] < 2 and processed_data.shape[1] >= 2:
                # If we have few rows but multiple columns, transpose the data
                # This handles cases where features are in columns instead of rows
                processed_data = processed_data.T
            
            # Very lenient validation - just ensure to have some data structure
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
        val1 = self.validate_species_name(species1)
        val2 = self.validate_species_name(species2)
        if not val1.get("valid"):
            return f"Invalid species 1: {val1.get('error','Unknown error')}"
        if not val2.get("valid"):
            return f"Invalid species 2: {val2.get('error','Unknown error')}"
        # Use possibly normalized names
        species1 = val1.get("species_name", species1)
        species2 = val2.get("species_name", species2)

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

    def _get_series_for_species(self, df, species):
        """Find species vector in rows or columns; return (values as float Series, axis_name)."""
        if species in df.index:
            s = df.loc[species].astype(float)
            return s, 'rows'
        if species in df.columns:
            s = df[species].astype(float)
            return s, 'cols'
        return None, None

    def _spearman_fast(self, x, y):
        """Spearman via rank+pearson (no scipy)."""
        xr = x.rank(method='average')
        yr = y.rank(method='average')
        return xr.corr(yr)

    def _binary_high_low(self, s):
        """Return boolean series of 'high' (>= median) and 'low' (< median)."""
        med = float(s.median())
        return (s >= med), (s < med), med

    def validate_species_name(self, species_name):
        """Validate and suggest corrections for species names"""
        if not species_name or not isinstance(species_name, str):
            return {"valid": False, "error": "Species name must be a non-empty string"}
    
        species_name = species_name.strip()
        if len(species_name) < 2:
            return {"valid": False, "error": "Species name too short"}
    
        # Check if species exists in current data
        if hasattr(self, 'current_analysis') and self.current_analysis:
            otu_data = self.current_analysis.get('otu_data')
            if otu_data is not None:
                # Check in both rows and columns
                available_species = list(otu_data.index) + list(otu_data.columns)
                if species_name not in available_species:
                    # Use existing suggestion method
                    suggestion = self._suggest_closest_taxa(
                        pd.DataFrame(index=available_species), 'index', species_name
                    )
                    return {
                        "valid": False, 
                        "error": f"Species '{species_name}' not found in dataset.{suggestion}",
                        "suggestions": suggestion
                    }
    
        return {"valid": True, "species_name": species_name}

    def infer_rule_evidence(self, species1, species2, high_thresh=0.7):
        """
        Build simple, human-readable rule evidence between two species from the OTU table.
        Returns dict with counts, correlations, rule_texts.
        """
        if not self.current_analysis or 'otu_data' not in self.current_analysis:
            return {'error': 'No OTU data available.'}

        df = self.current_analysis['otu_data']
        s1, ax1 = self._get_series_for_species(df, species1)
        s2, ax2 = self._get_series_for_species(df, species2)
        if s1 is None or s2 is None:
            # Try fuzzy suggestion if available
            hint = ''
            if hasattr(self, '_suggest_closest_taxa'):
                if s1 is None: hint += self._suggest_closest_taxa(df, df.index.name or 'index', species1)
                if s2 is None: hint += self._suggest_closest_taxa(df, df.index.name or 'index', species2)
            return {'error': f'Species not found in table. {hint}'.strip()}

        # align on common sample index
        s1, s2 = s1.align(s2, join='inner')
        if len(s1) < 3:
            return {'error': 'Not enough overlapping samples.'}

        # core stats
        spearman = float(self._spearman_fast(s1, s2))
        pearson  = float(s1.corr(s2))
        n = int(len(s1))

        # high/low evidence
        s1_hi, s1_lo, m1 = self._binary_high_low(s1)
        s2_hi, s2_lo, m2 = self._binary_high_low(s2)

        both_hi = int((s1_hi & s2_hi).sum())
        both_lo = int((s1_lo & s2_lo).sum())
        opp_hilo = int((s1_hi & s2_lo).sum())
        opp_lohi = int((s1_lo & s2_hi).sum())
        concordant = both_hi + both_lo
        discordant = opp_hilo + opp_lohi

        frac_conc = concordant / n
        frac_disc = discordant / n

        rule_texts = []
        if frac_conc >= high_thresh:
            rule_texts.append(
                f"If {species1} and {species2} are both high (≥ median) or both low in ≥{int(high_thresh*100)}% of samples "
                f"→ supports MUTUALISM/positive association (concordant={concordant}/{n})."
            )
        if frac_disc >= high_thresh:
            rule_texts.append(
                f"If {species1} is high when {species2} is low (or vice versa) in ≥{int(high_thresh*100)}% of samples "
                f"→ supports COMPETITION/negative association (discordant={discordant}/{n})."
            )
        if not rule_texts:

            rule_texts.append(
                f"Concordant={concordant}/{n} ({frac_conc:.2f}), Discordant={discordant}/{n} ({frac_disc:.2f}); "
                f"Spearman={spearman:.2f}, Pearson={pearson:.2f}."
            )

        return {
            'n_samples': n,
            'median': {species1: m1, species2: m2},
            'counts': {'both_high': both_hi, 'both_low': both_lo, 'hi_low': opp_hilo, 'low_hi': opp_lohi},
            'fractions': {'concordant': round(frac_conc, 3), 'discordant': round(frac_disc, 3)},
            'correlation': {'spearman': round(spearman, 3), 'pearson': round(pearson, 3)},
            'rule_texts': rule_texts
        }

    def proof_trace_for_interaction(self, species1, species2, inferred_type):
        """Produce a simple step-wise reasoning trace for UI/report."""
        ev = self.infer_rule_evidence(species1, species2)
        if 'error' in ev: return {'error': ev['error']}

        steps = []
        steps.append(f"1) Compute medians: {species1}={ev['median'][species1]:.3g}, {species2}={ev['median'][species2]:.3g}.")
        steps.append("2) Binarize each sample as High (≥ median) or Low (< median).")
        c = ev['counts']
        steps.append(f"3) Count patterns across {ev['n_samples']} samples: "
                 f"both-high={c['both_high']}, both-low={c['both_low']}, hi-low={c['hi_low']}, low-hi={c['low_hi']}.")
        steps.append(f"4) Concordant fraction={ev['fractions']['concordant']}, discordant={ev['fractions']['discordant']}; "
                 f"Spearman={ev['correlation']['spearman']}, Pearson={ev['correlation']['pearson']}.")
        if inferred_type.lower().startswith('mutual'):
            steps.append("5) Rule: high concordance → supports Mutualism/positive association.")
        elif inferred_type.lower().startswith('comp'):
            steps.append("5) Rule: high discordance → supports Competition/negative association.")
        else:
            steps.append("5) Rule: mixed patterns; label may stem from higher-order logic in R.")
        steps.append(f"6) Conclusion: supports '{inferred_type}' between {species1} and {species2}.")
        return {'steps': steps, 'evidence': ev}

    def network_null_significance(self, num_draws=200, random_state=42):
        """
        Show that inferred edges are not random by comparing mean |Spearman| on inferred pairs
        vs random pairs drawn from the same OTU table. Returns z-score and permutation p-value.
        """
        if not self.current_analysis or 'otu_data' not in self.current_analysis or 'interactions' not in self.current_analysis:
            return {'error': 'Need otu_data and interactions in memory.'}

        df = self.current_analysis['otu_data']
        inter = self.current_analysis['interactions']
        rng = np.random.default_rng(random_state)

        # Get species-axis names and vector access
        def get_vec(name):
            if name in df.index:  return df.loc[name].astype(float)
            if name in df.columns: return df[name].astype(float)
            return None

        # observed statistic: mean |Spearman| over inferred edges we actually have data for
        obs_vals = []
        for _, r in inter.iterrows():
            s1, s2 = str(r['sp1']), str(r['sp2'])
            v1, v2 = get_vec(s1), get_vec(s2)
            if v1 is None or v2 is None: continue
            v1, v2 = v1.align(v2, join='inner')
            if len(v1) >= 3:
                obs_vals.append(abs(self._spearman_fast(v1, v2)))
        if not obs_vals:
            return {'error': 'No overlapping species vectors for null test.'}
        obs = float(np.mean(obs_vals))

        # build random pairs from the same species pool
        species_pool = list(df.index if df.index.size >= df.columns.size else df.columns)
        species_pool = [str(s) for s in species_pool]
        num_pairs = len(obs_vals)

        null_stats = []
        for _ in range(num_draws):
            pairs = rng.choice(species_pool, size=(num_pairs, 2), replace=True)
            vals = []
            for a, b in pairs:
                va, vb = get_vec(a), get_vec(b)
                if va is None or vb is None: 
                    print(f"Warning: Missing data for species pair ({a}, {b})")
                    continue
                va, vb = va.align(vb, join='inner')
                if len(va) >= 3:
                    vals.append(abs(self._spearman_fast(va, vb)))
            if vals:
                null_stats.append(float(np.mean(vals)))

        if not null_stats:
            return {'error': 'Null distribution empty (species coverage issue).'}

        null_mean = float(np.mean(null_stats))
        null_std  = float(np.std(null_stats, ddof=1)) if len(null_stats) > 1 else 0.0
        z = (obs - null_mean) / (null_std + 1e-9)
        p = float((np.sum(np.array(null_stats) >= obs) + 1) / (len(null_stats) + 1))  # permutation p 

        return {
            'observed_mean_abs_spearman': round(obs, 3),
            'null_mean': round(null_mean, 3),
            'null_std': round(null_std, 3),
            'z_score': round(z, 3),
            'permutation_p_upper': round(p, 4),
            'num_inferred_pairs_tested': int(num_pairs),
            'num_null_draws': int(num_draws)
        }

    
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
        required_cols = {'sp1', 'sp2', 'lnk', 'comp'}
        missing = required_cols - set(map(str, interactions.columns))
        if missing:
            return {"error": f"Interactions frame missing columns: {', '.join(sorted(missing))}"}
        
        # Create NetworkX graph
        G = nx.from_pandas_edgelist(
            interactions, 
            source='sp1', 
            target='sp2', 
            edge_attr=['lnk', 'comp'],
            create_using=nx.DiGraph()
        )
        if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
            return {"error": "No nodes/edges to plot."}
        
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
        
        # Create edge traces by interaction type for legend
        edge_traces = []
        
        for interaction_type in [t for t in interactions['lnk'].dropna().unique()]:
            type_edges = interactions[interactions['lnk'] == interaction_type]
            type_x = []
            type_y = []
            hovertexts = []
            cds = []
            
            for _, row in type_edges.iterrows():
                s1, s2 = row['sp1'], row['sp2']
                if s1 not in pos or s2 not in pos:
                    continue
                x0, y0 = pos[s1]
                x1, y1 = pos[s2]
                type_x.extend([x0, x1, None])
                type_y.extend([y0, y1, None])
                comp_val = row['comp']
                try:
                    comp_txt = f"{float(comp_val):.3f}"
                except Exception:
                    comp_txt = str(comp_val)

                # one hover text per segment endpoint (Plotly expects per-point)
                ht = f"{s1} → {s2} ({interaction_type}, comp={comp_txt})"
                hovertexts.extend([ht, ht, None])

                # customdata per point (so clicks carry [sp1,sp2,lnk])
                cd = [s1, s2, interaction_type]
                cds.extend([cd, cd, None])

            if not type_x:
                continue
            
            edge_traces.append(
                go.Scatter(
                    x=type_x, y=type_y,
                    line=dict(width=3, color=interaction_colors.get(interaction_type, '#888888')),
                    hoverinfo='text',
                    mode='lines',
                    hovertext=hovertexts,
                    customdata=cds,
                    name=f"{interaction_type} ({len(type_edges)})",
                    showlegend=True
                ))
            
        
        # This creates node trace with better styling
        node_x = []
        node_y = []
        node_text = []
        node_info = []
        node_sizes = []
        
        for node in G.nodes():
            if node not in pos:
                continue
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            
            # Use actual species name (truncate if too long)
            label = str(node)
            node_text.append(label if len(label) <= 15 else label[:12] + "…")

            indeg = G.in_degree(node)
            outdeg = G.out_degree(node)
            total = indeg + outdeg
            node_sizes.append(max(20, min(50, 20 + total * 3)))

            node_info.append(
                f"<b>Species: {node}</b><br>"
                f"Incoming: {indeg}<br>Outgoing: {outdeg}<br>Total degree: {total}"
            )
        
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
        
        # To Create figure with better layout
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
                               text="Hover over nodes for details and explanation of interaction.",
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
        # Quick history access via chat (natural language)
        
        # Dataset comparison via chat (e.g., 'compare interactions_enhanced vs interactions_short')
        if 'compare' in query_lower:
            def _pick_two(hintA, hintB):
                a_df, a_name = self._load_csv_from_known_locations(hintA)
                b_df, b_name = self._load_csv_from_known_locations(hintB)
                return a_df, a_name, b_df, b_name

            if 'interaction' in query_lower:
                # parse specific names if provided around 'vs'
                a_hint, b_hint = 'interactions_enhanced', 'interactions_short'
                if ' vs ' in query_lower:
                    parts = [p.strip() for p in query_lower.split(' vs ')]
                    if len(parts) == 2:
                        # try to extract filename-like tokens
                        a_hint = parts[0].replace('compare', '').strip() or a_hint
                        b_hint = parts[1].strip() or b_hint
                a_df, a_name, b_df, b_name = _pick_two(a_hint, b_hint)
                if a_df is None or b_df is None:
                    return f"Couldn't find both interaction tables ('{a_hint}' and '{b_hint}'). Try uploading them first."
                res = self.compare_interactions(a_df, b_df)
                txt = [
                    f"📊 Interaction tables compared: **{a_name}** vs **{b_name}**",
                    f"- Edges: {res['a_edges']} vs {res['b_edges']}",
                    f"- Overlap: {res['overlap_edges']} | Jaccard={res['jaccard']:.3f}",
                ]
                if res['added_examples']:
                    ex = ', '.join([f"{x[0]}—{x[1]}" for x in res['added_examples'][:5]])
                    txt.append(f"- Examples only in {b_name}: {ex} ...")
                if res['removed_examples']:
                    ex = ', '.join([f"{x[0]}—{x[1]}" for x in res['removed_examples'][:5]])
                    txt.append(f"- Examples only in {a_name}: {ex} ...")
                if res['label_changes']:
                    k = next(iter(res['label_changes'].keys()))
                    ch = res['label_changes'][k]
                    txt.append(f"- Label change example {k[0]}—{k[1]}: {ch['a']} → {ch['b']}")
                return "\n".join(txt)

            if 'taxonom' in query_lower:
                a_df, a_name = self._load_csv_from_known_locations('taxonomy')
                b_df, b_name = self._load_csv_from_known_locations('taxonomy_table')
                if a_df is None or b_df is None:
                    return "Couldn't find two taxonomy files to compare. Upload or place them in the data folder."
                res = self.compare_taxonomy(a_df, b_df)
                if 'error' in res:
                    return res['error']
                lines = [f"🧬 Taxonomy compare: **{a_name}** vs **{b_name}** (columns: {', '.join(res['columns_compared'])})"]
                for lvl in res['levels']:
                    lines.append(f"- {lvl['relation']}: pairs {lvl['a_pairs']} vs {lvl['b_pairs']}, overlap {lvl['overlap']}, Jaccard={lvl['jaccard']:.3f}")
                return "\n".join(lines)

            if 'mapping' in query_lower or 'species map' in query_lower:
                a_df, a_name = self._load_csv_from_known_locations('mapping')
                b_df, b_name = self._load_csv_from_known_locations('augmented')
                if a_df is None or b_df is None:
                    return "Couldn't find two mapping files to compare. Upload mapping and augmented mapping first."
                res = self.compare_mappings(a_df, b_df)
                if 'error' in res:
                    return res['error']
                return (f"🧾 Mapping compare: **{a_name}** vs **{b_name}**\n"
                        f"- Keys {res['a_keys']} vs {res['b_keys']}, overlap {res['overlap']}, Jaccard={res['jaccard']:.3f}")
        if any(term in query_lower for term in ['history', 'my history', 'queries this session', 'show history', 'show my history']):
            if not self.chat_history:
                return "No history yet. Ask me something!"
            lines = []
            for i, h in enumerate(self.chat_history[-10:], 1):
                u = h.get('user','').strip()
                b = h.get('bot','').strip()
                ts = h.get('timestamp','')
                lines.append(f"{i}. [{ts}] You: {u}\n   Bot: {b[:180]}{'...' if len(b)>180 else ''}")
            return "**Recent queries (last 10):**\n\n" + "\n\n".join(lines)
    
        
        # Handle shutdown request
        if query_lower in ['shutdown', 'exit', 'quit']:
            return "SHUTDOWN_REQUESTED"
        
        # Handle static Q&A first ( which work regardless of uploaded files)
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

        if 'help' in query_lower:
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
- "How does InfIntE work?"
- "Show me interaction types"
- "List all species in the data"
- "What families are present?"
            """

        if any(term in query_lower for term in ['delete selected csv', 'delete selected csvs', 'remove selected csv', 'remove selected csvs']):
            return "DELETE_SELECTED_CSVS"
        # Handle CSV data queries only after static Q&A
        csv_response = self.handle_csv_data_query(query_lower)
        if csv_response:
            return csv_response
        
        # Analysis Queries 
        # Dataset filtering commands 
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
                    
                    if not self.current_analysis:
                        return "No analysis data available. Please upload and analyze data first."
                    if 'interactions' not in self.current_analysis:
                        return "No interaction data found in current analysis."

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

            except KeyError as e:
                if 'current_analysis' in str(e):
                    return "No analysis data available. Please upload and analyze data first."
                elif 'comp' in str(e):
                    return "Compression data not available for abundance filtering. The analysis may not have completed successfully."
                else:
                    return f"Missing data field: {str(e)}"
            except AttributeError:
                return "No analysis data available. Please upload and analyze data first."

            except Exception as e:
                return f"Error processing filter request: {str(e)}. Please try uploading and analyzing data first."
        
        # ---- PDF Search Queries ----
        explicit_pdf_terms = ["pdf","according to the paper", "paper", "from the pdf", "paper says", "document states", "literature shows", "study mentions", "research paper", "in pdf"]
        if any(term in query_lower for term in explicit_pdf_terms) and self.uploaded_files['pdf_files']:
            return self.search_across_all_pdfs(query)
        
        # ---- CSV Data Queries ----
        # Handles questions about uploaded CSV data
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
        pdf_hints = [
        "pdf", "paper", "document", "literature", "study",
        "according to the pdf", "according to pdf",
        "according to the paper", "from the pdf", "from pdf",
        "paper says", "document states", "literature shows",
        "study mentions", "in the pdf", "in pdf"
        ]
        if self.uploaded_files['pdf_files'] and any(term in query_lower for term in pdf_hints):
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
            return 'otu_table'  
        
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
            
            return None  # If no match found
            
        except Exception as e:
            return f"Error processing CSV query: {str(e)}"
    
    def handle_intelligent_taxonomy_query(self, df: pd.DataFrame, filename: str, query_lower: str) -> str:
        """Handle taxonomy queries with intelligent parsing and filtering"""
        import re
        
        dna_hits = re.findall(r'[ACGTNacgtn]{20,}', query_lower.replace(" ", ""))
        if dna_hits:
            seq = max(dna_hits, key=len)  # take the longest stretch if multiple present
            return self.lookup_full_taxonomy_for_otu_sequence(seq)
            
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
        #   A) "<child> of <parent> <name>"                e.g., "phylum of order Micrococcales"
        #   B) "<parent> <name> has which <child>"         e.g., "order Micrococcales has which phylum"
        #   C) "which <child> are in <parent> <name>"      e.g., "which families are in order Micrococcales"
        #   D) "<child> in <parent> <name>"                e.g., "families in order Micrococcales"
        #   E) "show <child> for <parent> <name>"          e.g., "show families for order Micrococcales"
        rel_patterns = [
            (r'\b(kingdom|phylum|class|order|family|genus)\s+of\s+(kingdom|phylum|class|order|family|genus)\s+([\w\-]+)\b', 'A'),
            (r'\b(kingdom|phylum|class|order|family|genus)\s+([\w\-]+)\s+has\s+which\s+(kingdom|phylum|class|order|family|genus)\b', 'B'),
            (r'\bwhich\s+(kingdom|phylum|class|order|family|genus)\s+are\s+in\s+(kingdom|phylum|class|order|family|genus)\s+([\w\-]+)\b', 'C'),
            (r'\b(kingdom|phylum|class|order|family|genus)\s+in\s+(kingdom|phylum|class|order|family|genus)\s+([\w\-]+)\b', 'D'),
            (r'\bshow\s+(kingdom|phylum|class|order|family|genus)\s+for\s+(kingdom|phylum|class|order|family|genus)\s+([\w\-]+)\b', 'E'),
        ]
        for pat, kind in rel_patterns:
            _m = re.search(pat, query_lower)
            if not _m:
                continue
            if kind == 'A':
                child_level  = _m.group(1).lower()
                parent_level = _m.group(2).lower()
                parent_name  = _m.group(3).strip()
            elif kind == 'B':
                parent_level = _m.group(1).lower()
                parent_name  = _m.group(2).strip()
                child_level  = _m.group(3).lower()
            else:  
                child_level  = _m.group(1).lower()
                parent_level = _m.group(2).lower()
                parent_name  = _m.group(3).strip()

            if parent_level not in taxonomy_levels or child_level not in taxonomy_levels:
                return f"Your taxonomy file doesn’t have '{child_level}' and '{parent_level}' columns."

            parent_col = taxonomy_levels[parent_level]
            child_col  = taxonomy_levels[child_level]

            # Exact match first, then partial
            rows = df[df[parent_col].astype(str).str.lower() == parent_name.lower()]
            if rows.empty:
                rows = df[df[parent_col].astype(str).str.contains(parent_name, case=False, na=False)]

            if not rows.empty:
                vals = sorted(rows[child_col].dropna().astype(str).unique())
                if len(vals) == 1:
                    return f"**{child_col} of {parent_col} '{parent_name}':**\n\n• {vals[0]}"
                bullets = "\n".join(f"• {v}" for v in vals[:20])
                extra = "" if len(vals) <= 20 else f"\n… and {len(vals)-20} more"
                return f"**{child_col} of {parent_col} '{parent_name}':**\n\n{bullets}{extra}"

            # Nothing found -> suggest closest parent names if helper exists
            try:
                hint = self._suggest_closest_taxa(df, parent_col, parent_name)
            except Exception:
                hint = ""
            return f"No records found for {parent_level} '{parent_name}' in {filename}.{hint}"
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
                                return f"No records found for {child_level} '{child_name}' in {filename}." + self._suggest_closest_taxa(df, child_col, child_name)
        
        # Handle simple listing queries
        simple_queries = {
            'phyla': 'phylum', 'phylum': 'phylum', 'what phyla': 'phylum',
            'families': 'family', 'family': 'family', 'what families': 'family', 'show families': 'family',
            'classes': 'class', 'class': 'class', 'what classes': 'class',
            'orders': 'order', 'order': 'order', 'what orders': 'order',
            'genus': 'genus', 'genera': 'genus', 'what genus': 'genus', 'list genus': 'genus', 'show genus': 'genus',
            'species': 'species', 'what species': 'species', 'list species': 'species',
            'kingdom': 'kingdom', 'kingdoms': 'kingdom', 'what kingdoms': 'kingdom', 'show kingdoms': 'kingdom',
        }
        
        for query_term, level in simple_queries.items():
            if query_term in query_lower:
                if level in taxonomy_levels:
                    return self.handle_taxonomy_query(df, filename, level)
        
        # Handle taxonomy summary
        if any(term in query_lower for term in ['taxonomy', 'taxonomic', 'classification']):
            return self.handle_taxonomy_summary(df, filename)
        
        # Default fallback
        return f"I couldn't understand your taxonomy query. Try asking (also works: 'phylum of order <Name>'):\n" + \
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
        
        for i, result in enumerate(results[:5], 1):  
            response += f"**{i}. From {result['filename']}:**\n"
            response += f"{result['context']}\n"
            response += f"*Relevance: {result['relevance_score']:.2f}*\n\n"
        
        if len(results) > 5:
            response += f"*...and {len(results) - 5} more results. Try a more specific query for better results.*"
        
        response = response.strip()
        return response

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
            
            # quick path: if user pasted a DNA-like sequence and asks for taxonomy/all ranks
            dna_like = re.search(r'\b[ACGTNacgtn]{20,}\b', query_lower, flags=re.IGNORECASE)
            asks_tax = any(w in query_lower for w in [
                'taxonomy', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'
            ])
            if dna_like and asks_tax:
                seq = dna_like.group(0)
                return self.lookup_full_taxonomy_for_otu_sequence(seq)

            dna_prompt = re.search(r'species\s+([ACGTNacgtn]{20,})', query_lower, flags=re.IGNORECASE)
            if dna_prompt and asks_tax:
                return self.lookup_full_taxonomy_for_otu_sequence(dna_prompt.group(1))

            # Try to answer any question about the CSV data intelligently
            return self.handle_intelligent_data_query(df, filename, csv_type, query_lower)
            
        except Exception as e:
            return f"Error processing query: {str(e)}"
    
    def handle_otu_id_lookup(self, query_lower: str) -> str:
        """Handle OTU ID to species name lookup queries and species to OTU ID queries"""
        import re

        # Check for species name to OTU ID queries 
        species_patterns = [
            r'otu.{0,5}id.{0,5}for.{0,5}(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'what.{0,5}is.{0,5}the.{0,5}otu.{0,5}id.{0,5}for.{0,5}(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'otu.{0,5}id.{0,5}of.{0,5}(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'find.{0,5}otu.{0,5}id.{0,5}of.{0,5}(?:the.{0,5})?(?:species.{0,5}name.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)',
            r'find.{0,5}otu.{0,5}id.{0,5}(?:for.{0,5})?([A-Za-z][A-Za-z0-9_\-]+)(?:\?|$)'
        ]
        
        for pattern in species_patterns:
            match = re.search(pattern, query_lower)
            if match:
                species_name = match.group(1).strip()
                return self.lookup_species_to_otu(species_name)
        
        
        # Check for numeric ID queries 
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
        
        # Check for OTU sequence queries 
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
            r'(?:(?:kingdom|phylum|class|order|family|genus|species)\s+)?([A-Za-z][A-Za-z0-9_\-]+)\s+(?:has|belongs?\s+to|is\s+in)\s+which\s+(?:kingdom|phylum|class|order|family|genus)',
            r'genus\s+([A-Za-z][A-Za-z0-9_\-]+)\s+(?:has|belongs?\s+to|is\s+in)\s+which\s+(?:phylum|class|order|family|kingdom)',
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
        END = r'(?:[?.!]|$)'
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

        unique_name_by_otu_patterns = [
            r'(?:what.{0,5}is.{0,5})?unique.{0,5}name.{0,5}(?:for|of).{0,5}otu.{0,5}id.{0,5}`?([A-Za-z0-9_\-]+)`?(?:\?|$)',
            r'unique.{0,5}name.{0,5}for.{0,5}otu.{0,5}id.{0,5}`?([A-Za-z0-9_\-]+)`?(?:\?|$)',
            r'unique.{0,5}name.{0,5}of.{0,5}otu.{0,5}id.{0,5}`?([A-Za-z0-9_\-]+)`?(?:\?|$)',
        ]

        for pattern in unique_name_by_otu_patterns:
            match = re.search(pattern, query_lower)
            if match:
                otu_id = match.group(1).strip()
                return self.lookup_unique_name_by_otu_id(otu_id)

        
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
        
        otu_id = otu_matches[0]  # THis takes the first/longest match
        
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
    
    def lookup_unique_name_by_otu_id(self, otu_id: str) -> str:
        """Look up unique name for a given OTU ID"""
        import re

        # Normalize the input (strip whitespace/backticks)
        otu_id_norm = str(otu_id).strip().strip('`')

        # Search through uploaded CSV files
        for filename, file_info in self.uploaded_files['csv_files'].items():
            df = file_info['data']

        # Need at least otu_id + species_name; unique_name is optional
            if 'otu_id' in df.columns and 'species_name' in df.columns:
                otu_col = df['otu_id'].astype(str)

            # Exact match (case-insensitive)
                exact = df[otu_col.str.casefold() == otu_id_norm.casefold()]

                if len(exact) > 0:
                    row = exact.iloc[0]
                    species = row.get('species_name', '')
                    unique = row.get('unique_name') if 'unique_name' in df.columns else None

                    return (
                        "**Unique Name Lookup Result:**\n\n"
                        f"🔬 **Species Name:** **{species}**\n"
                        f"🧬 **OTU ID:** `{otu_id_norm}`\n"
                        + (f"🏷️ **Unique Name:** **{str(unique)}**\n\n" if unique is not None and str(unique).strip() != "" else "\n")
                        + f"*Found in {filename}*"
                    )

            # Partial match (helpful if OTU IDs are long strings)
                partial = df[otu_col.str.casefold().str.startswith(otu_id_norm.casefold())]

                if len(partial) > 0:
                    results = []
                    for _, row in partial.head(3).iterrows():  # Show top 3 matches
                        pid = str(row.get('otu_id', ''))
                        pspp = row.get('species_name', '')
                        punq = row.get('unique_name', '') if 'unique_name' in df.columns else ''
                        if punq and str(punq).strip() != "":
                            results.append(f"• **{pspp}** — OTU `{pid}` → {punq}")
                        else:
                            results.append(f"• **{pspp}** — OTU `{pid}`")

                    return (
                        "**Partial OTU ID Matches Found:**\n\n"
                        f"🔍 **Search term:** `{otu_id_norm}`\n\n"
                        + "\n".join(results)
                        + f"\n\n*Found in {filename}*"
                    )

        return (
            "**OTU ID Not Found:**\n\n"
            f"🧬 **Query OTU ID:** `{otu_id_norm}`\n"
            "❌ No matching record with a unique name found in uploaded files.\n\n"
            "Please ensure your CSV has 'otu_id', 'species_name', and (optionally) 'unique_name' columns."
        )

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
                           f"{additional_info}\n" + \
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
                           f"{additional_info}\n" + \
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
        """Look up taxonomy information for any rank pair, plus legacy fallbacks."""
        import re

        q = (query_lower or "").strip().lower()
        token = (search_term or "").strip().lower()

       
        m = re.search(
            r'\b(?:(kingdom|phylum|class|order|family|genus|species)\s+)?'   
            r'([a-z0-9_\-]+)\s+'                                             
            r'(?:has|belongs?\s+to|is\s+in)\s+which\s+'                      
            r'(kingdom|phylum|class|order|family|genus)\b',                 
            q
        )
        if m:
            src_rank, name, dst_rank = m.groups()
            return self._dispatch_taxonomy_pair(name, src_rank, dst_rank)

        # B) "what is the <rank> of <rank?> <name>"
        m = re.search(
            r'\b(?:what\s+is\s+the\s+)?'
            r'(kingdom|phylum|class|order|family|genus)\s+of\s+'
            r'(?:(kingdom|phylum|class|order|family|genus)\s+)?'
            r'([a-z0-9_\-]+)\b',
            q
        )
        if m:
            dst_rank, src_rank, name = m.groups()
            return self._dispatch_taxonomy_pair(name, src_rank, dst_rank)

        # --- Legacy specific fallbacks (This is kept for backward compatibility) ---
        if 'family name of genus' in q:
            return self.lookup_family_of_genus(token)
        elif 'genus of' in q and 'family' in q:
            return self.lookup_genus_of_family(token)
        elif 'phylum of order' in q:
            return self.lookup_phylum_of_order(token)
        elif 'phylum of genus' in q:
            return self.lookup_phylum_of_genus(token)
        elif 'class of genus' in q:
            return self.lookup_class_of_genus(token)
        elif 'order of genus' in q:
            return self.lookup_order_of_genus(token)
        elif 'order of phylum' in q:
            return self.lookup_order_of_phylum(token)
        elif 'kingdom of' in q:
            return self.lookup_kingdom_of_taxon(token)
        elif 'species of genus' in q:
            return self.lookup_species_of_genus(token)
        elif 'species' in q and any(r in q for r in ['order', 'family', 'class', 'phylum', 'kingdom']):
            return self.lookup_species_by_taxonomy(token, q)

        # Default
        return self.lookup_family_of_genus(token)

    def _dispatch_taxonomy_pair(self, name, src_rank, dst_rank) -> str:
        """Route (source_rank -> target_rank) to existing helpers; else fallback to generic table lookup."""
        name = (name or "").strip()
        src = (src_rank or "").lower() or None
        dst = (dst_rank or "").lower()

        # common fast path
        if dst == "kingdom":
            return self.lookup_kingdom_of_taxon(name)

        # use existing specific helpers when available
        mapping = {
        ("genus",  "family"):  self.lookup_family_of_genus,
        ("genus",  "order"):   self.lookup_order_of_genus,
        ("genus",  "class"):   self.lookup_class_of_genus,
        ("genus",  "phylum"):  self.lookup_phylum_of_genus,
        ("family", "genus"):   self.lookup_genus_of_family,
        ("order",  "phylum"):  self.lookup_phylum_of_order,
        ("phylum", "order"):   self.lookup_order_of_phylum,
    }

        if src and (src, dst) in mapping:
            return mapping[(src, dst)](name)

        if not src:
            for try_src in ("genus", "family", "order", "class", "phylum"):
                if (try_src, dst) in mapping:
                    try:
                        return mapping[(try_src, dst)](name)
                    except Exception:
                        pass  # try next

        # generic table-driven answer (works for any pair found in taxonomy table)
        return self._lookup_rank_generic(name, dst, src)

    def _lookup_rank_generic(self, token, target_rank, source_rank=None) -> str:
        """Generic, table-driven answer for any (source -> target) pair using uploaded taxonomy CSV."""
        token_l = (token or "").strip().lower()
        target = (target_rank or "").lower()
        source = (source_rank or "").lower() if source_rank else None

        # To find a taxonomy dataframe among uploaded CSVs
        csvs = (self.uploaded_files or {}).get("csv_files", {})
        best_df, best_name, cols = None, None, None
        for fname, meta in csvs.items():
            df = meta.get("data") if isinstance(meta, dict) else None
            if df is None or getattr(df, "empty", False):
                continue
            cols_lower = {c.lower(): c for c in df.columns}
            if len({"kingdom","phylum","class","order","family","genus","species"} & set(cols_lower.keys())) >= 3:
                best_df, best_name, cols = df, fname, cols_lower
                break
        if best_df is None:
            return "No taxonomy table found. Please upload a taxonomy CSV with columns like Genus/Family/Order/Phylum."

        candidates = best_df
        if source and source in cols:
            candidates = candidates[candidates[cols[source]].astype(str).str.strip().str.lower() == token_l]
        else:
            hit = False
            for r in ("genus","species","family","order","class","phylum","kingdom"):
                if r in cols:
                    mask = candidates[cols[r]].astype(str).str.strip().str.lower() == token_l
                    if mask.any():
                        candidates = candidates[mask]
                        source = r
                        hit = True
                    break
        if not hit:
            return f"Could not find '{token}' in the taxonomy table."

        if target not in cols:
            return f"'{target_rank}' column not present in the taxonomy table."

        values = candidates[cols[target]].dropna().astype(str).str.strip().unique().tolist()
        if not values:
            return f"No {target_rank} found for {source or 'taxon'} '{token}'."

        pretty_src = (source or "taxon").capitalize()
        pretty_dst = target_rank.capitalize()
        src_val    = token
        dst_val    = values[0] if len(values) == 1 else ", ".join(values[:10]) + (" …" if len(values) > 10 else "")

        return (
            f"<div><strong>Taxonomy Lookup Result:</strong><br>"
            f"🧬 <strong>{pretty_src}:</strong> {src_val}<br>"
            f"🌿 <strong>{pretty_dst}:</strong> {dst_val}<br>"
            f"<em>*Found in {best_name}*</em></div>"
        )

    
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
        # First, it finds the species name for this OTU sequence
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
        
        # Searches taxonomy files for this species
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
    

    def _is_dna_series(self, s, min_ratio=0.7, min_len=20):
        import re
        vals = s.dropna().astype(str).head(200)
        if vals.empty:
            return False
        dna = vals.str.fullmatch(r'[ACGTNacgtn]+', na=False)
        long = vals.str.len().ge(min_len)
        return (dna & long).mean() >= min_ratio

    def _normcolmap(self, df):
        """Case-insensitive map: lower -> original col name"""
        return {c.lower(): c for c in df.columns}

    def _looks_like_sequence_col(self, colname: str) -> bool:
        colname_l = colname.lower()
        
        return any(k in colname_l for k in ["sequence", "dna", "otu", "id"])

    def _first_taxonomy_row_for_species(self, species: str):
        """Search any uploaded CSV that looks like a taxonomy table for a row about this species."""
        wanted = ["kingdom","phylum","class","order","family","genus","species"]
        for fname, finfo in (self.uploaded_files.get("csv_files") or {}).items():
            df = finfo.get("data")
            if df is None or not isinstance(df, pd.DataFrame):
                continue
            cmap = self._normcolmap(df)
            have = [r for r in wanted if r in cmap]
            if len(have) < 3:
                continue  

            # prefer matching by Species if present
            if "species" in cmap:
                rows = df[df[cmap["species"]].astype(str).str.lower().str.contains(species.lower(), na=False)]
                if not rows.empty:
                    return rows.iloc[0], fname

            # otherwise: try any rank columns
            mask = pd.Series(False, index=df.index)
            for r in have:
                mask |= df[cmap[r]].astype(str).str.lower().str.contains(species.lower(), na=False)
            rows = df[mask]
            if not rows.empty:
                return rows.iloc[0], fname
        return None, None

    def _species_from_sequence(self, seq: str):
        """
        Try to resolve a species name from a sequence by:
        1) exact/loose match against obvious sequence/OTU columns
        2) using an explicit mapping table if present (columns: otu_id -> species_name)
        Returns (species_name or None, provenance_info)
        """
        seq_u = seq.upper()

        # 1) Search any 'sequence-like' column
        for fname, finfo in (self.uploaded_files.get("csv_files") or {}).items():
            df = finfo.get("data")
            if df is None or not isinstance(df, pd.DataFrame):
                continue
            cmap = self._normcolmap(df)
            for col in df.columns:
                is_seq_col = self._looks_like_sequence_col(col) or self._is_dna_series(df[col])
                if not is_seq_col:
                    continue
                s = df[col].astype(str)
                
                # exact match first
                su = s.str.upper()
                hits = df[su == seq_u]
                if not hits.empty:
                    # if a species column is here, take it; or fall back to best of taxonomic columns
                    
                    if "species" in cmap:
                        return str(hits.iloc[0][cmap["species"]]).strip(), f"{fname}:{col}"
                    # try taking the most specific of genus/species-like columns
                    for k in ["genus","family","order","class","phylum","kingdom"]:
                        if k in cmap:
                            return str(hits.iloc[0][cmap[k]]).strip(), f"{fname}:{col}"
                
                hits = df[su.str.contains(seq_u[:50], na=False, regex=False)]
                if not hits.empty:
                    if "species" in cmap:
                        return str(hits.iloc[0][cmap["species"]]).strip(), f"{fname}:{col}~contains"
                    for k in ["genus","family","order","class","phylum","kingdom"]:
                        if k in cmap:
                            return str(hits.iloc[0][cmap[k]]).strip(), f"{fname}:{col}~contains"

        # 2) Tries explicit mapping tables (otu_id -> species_name)
        for fname, finfo in (self.uploaded_files.get("csv_files") or {}).items():
            df = finfo.get("data")
            if df is None or not isinstance(df, pd.DataFrame):
                continue
            cmap = self._normcolmap(df)
            if "otu_id" in cmap and "species_name" in cmap:
                s = df[cmap["otu_id"]].astype(str).str.upper()
                hits = df[s == seq_u]
                if not hits.empty:
                    return str(hits.iloc[0][cmap["species_name"]]).strip(), f"{fname}:mapping"
                hits = df[s.str.contains(seq_u[:50], na=False, regex=False)]
                if not hits.empty:
                    return str(hits.iloc[0][cmap["species_name"]]).strip(), f"{fname}:mapping~contains"

        return None, None

    def lookup_full_taxonomy_for_otu_sequence(self, sequence: str) -> str:
        """
        Main entry: DNA/OTU sequence -> (Species) -> full ranks.
        Works with taxonomy_table_short.csv (Kingdom..Genus..Species) or any similar file.
        """
        if not sequence or len(sequence.strip()) < 10:
            return "Please provide a valid DNA/OTU sequence."

        seq = re.sub(r"[^ACGTNacgtn]", "", sequence).upper()
        if len(seq) < 10:
            return "That doesn't look like a valid DNA-like sequence."

        # find species (or closest label) from uploaded CSVs
        species, prov = self._species_from_sequence(seq)
        if not species:
            return (
                "**OTU/Sequence not found in uploaded data.**\n\n"
                f"🧬 Query: `{seq[:60]}…`\n"
                "Please upload a mapping (e.g., columns **otu_id, species_name**) or a CSV that contains both the sequence and taxonomy."
            )

        # get taxonomy row
        row, src = self._first_taxonomy_row_for_species(species)
        if row is None:
            return (
                f"**Species resolved:** {species}\n"
                "But I couldn't find a taxonomy row (Kingdom..Genus) in your uploaded CSVs. "
                "Upload a taxonomy CSV (e.g., `taxonomy_table_short.csv`) with columns like "
                "`Kingdom, Phylum, Class, Order, Family, Genus, Species`."
            )

        def getv(keys):
            for k in keys:
                if k in row.index:
                    v = str(row[k]).strip()
                    if v and v != 'nan':
                        return v
            return "—"

        K = getv(["Kingdom","kingdom"])
        P = getv(["Phylum","phylum"])
        C = getv(["Class","class"])
        O = getv(["Order","order"])
        F = getv(["Family","family"])
        G = getv(["Genus","genus"])
        S = getv(["Species","species"])

        out = [
            "**Full Taxonomy**",
            f"🧬 **Sequence:** `{seq[:60]}…`",
            f"🔬 **Species:** {S if S != '—' else species}",
            "",
            f"• **Kingdom:** {K}",
            f"• **Phylum:** {P}",
            f"• **Class:** {C}",
            f"• **Order:** {O}",
            f"• **Family:** {F}",
            f"• **Genus:** {G}",    
        ]
        src_hint = f"\n*from {src}*" if src else ""
        if prov:
            out.append(f"\n*(resolved via {prov})*")
        if src_hint:
            out.append(src_hint)
        return "\n".join(out)

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
        
        # Extracts column names and data types
        columns = df.columns.tolist()
        numeric_columns = df.select_dtypes(include=['number']).columns.tolist()
        text_columns = df.select_dtypes(include=['object']).columns.tolist()
        
        # 1. Handles specific column queries (but avoid generic column names in species queries)
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
        
        for col in text_columns[:5]:  # Limits to first 5 columns
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
        
        for col in numeric_columns[:5]:  # Limits to first 5 numeric columns
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
            return (f"**Search Help for {filename}:**\n" "Please put your search term inside quotes.\n" "Example: 'find \"microvirga\"' or 'search \"bacteria\"'")
        
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

# This Initializes chatbot
chatbot = MicrobiomeXAIChatbot()

@app.route('/')
def index():
    """Main chatbot interface"""
    return render_template('chatbot.html')

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

    if response == "DELETE_SELECTED_CSVS":
        return jsonify({
            'response': 'Deleting selected CSV files...',
            'delete_csvs': True
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
        print("Analysis started...")
        
        result = chatbot.analyze_microbiome_data(otu_data, depth)
        print("Preprocessing complete...")
        print("R analysis running...")
        print("Results processed...")
        result['file_stored'] = True
        return jsonify(result)
               
    except Exception as e:
        return jsonify({'error': str(e)}), 400

@app.route('/network')
def network():
    """Generate network visualization"""
    result = chatbot.generate_network_visualization()
    return jsonify(result)


@app.route('/analyze_multiple', methods=['POST'])
def analyze_multiple():
    # Analyze multiple CSVs sent as JSON payload 
    # Stores all CSVs and returns a summary message. 
    try:
        payload = request.json or {}
        files = payload.get('files', [])
        if not files or not isinstance(files, list):
            return jsonify({'error': 'No files provided'}), 400

        ingested = []
        type_counts = {}
        for item in files:
            fname = item.get('filename', 'uploaded.csv')
            df = pd.DataFrame(item.get('otu_data', {}))
            chatbot.uploaded_files['csv_files'][fname] = {
                'path': fname,
                'data': df,
                'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            }
            try:
                csv_type = chatbot.detect_csv_type(df, fname)
            except Exception:
                csv_type = 'generic'
            ingested.append({'filename': fname, 'rows': len(df), 'cols': len(df.columns), 'type': csv_type})
            type_counts[csv_type] = type_counts.get(csv_type, 0) + 1

        msg = "✅ Ingested {} file(s). ".format(len(ingested))
        if type_counts:
            msg += " Types: " + ", ".join(f"{k}×{v}" for k, v in type_counts.items()) + "."
        msg += " You can now ask for summaries, filters, or visualize the network (if interaction data)."

        return jsonify({'success': True, 'message': msg, 'ingested': ingested})
    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/analyze_pairwise', methods=['POST'])
def analyze_pairwise():
    # Pairwise analysis for exactly two CSVs. If types match (interactions/taxonomy/mapping),
    # use existing comparators; otherwise just ingest both.
    try:
        payload = request.json or {}
        files = payload.get('files', [])
        if not files or len(files) != 2:
            return jsonify({'error': 'Exactly two files are required for pairwise analysis'}), 400

        dfs = []
        meta = []
        for item in files:
            fname = item.get('filename', 'uploaded.csv')
            df = pd.DataFrame(item.get('otu_data', {}))
            chatbot.uploaded_files['csv_files'][fname] = {
                'path': fname,
                'data': df,
                'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            }
            try:
                csv_type = chatbot.detect_csv_type(df, fname)
            except Exception:
                csv_type = 'generic'
            dfs.append(df)
            meta.append({'filename': fname, 'type': csv_type, 'rows': len(df), 'cols': len(df.columns)})

        t1, t2 = meta[0]['type'], meta[1]['type']

        if t1 == t2 == 'interactions':
            res = chatbot.compare_interactions(dfs[0], dfs[1])
            return jsonify({'success': True, 'pairwise_type': 'interactions', 'message': 'Compared interaction tables.', 'comparison': res, 'files': meta})
        elif t1 == t2 == 'taxonomy':
            res = chatbot.compare_taxonomy(dfs[0], dfs[1])
            return jsonify({'success': True, 'pairwise_type': 'taxonomy', 'message': 'Compared taxonomy tables.', 'comparison': res, 'files': meta})
        elif t1 == t2 == 'mapping':
            res = chatbot.compare_mappings(dfs[0], dfs[1])
            return jsonify({'success': True, 'pairwise_type': 'mapping', 'message': 'Compared mapping tables.', 'comparison': res, 'files': meta})
        else:
            return jsonify({
                'success': True,
                'pairwise_type': 'ingested',
                'message': 'Files ingested. Types differ or are not pairwise-comparable. You can still query each or run single-file analyses.',
                'files': meta
            })
    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/preview_csvs', methods=['POST'])
def preview_csvs():
    # Return a preview (first N rows) for selected CSVs already uploaded.
    try:
        payload = request.json or {}
        names = payload.get('filenames', [])
        if not names:
            return jsonify({'success': False, 'error': 'No filenames provided'}), 400

        previews = []
        for name in names:
            info = chatbot.uploaded_files['csv_files'].get(name)
            if not info:
                match = None
                for k in chatbot.uploaded_files['csv_files'].keys():
                    if k.lower() == str(name).lower():
                        match = chatbot.uploaded_files['csv_files'][k]
                        name = k
                        break
                info = match
            if not info:
                previews.append({'filename': name, 'error': 'Not found'})
                continue
            df = info.get('data')
            try:
                csv_type = chatbot.detect_csv_type(df, name)
            except Exception:
                csv_type = 'generic'
            head = df.head(30) if hasattr(df, 'head') else df
            previews.append({
                'filename': name,
                'csv_type': csv_type,
                'columns': list(head.columns) if hasattr(head, 'columns') else [],
                'rows': head.to_dict(orient='records') if hasattr(head, 'to_dict') else []
            })

        return jsonify({'success': True, 'previews': previews})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 400

@app.route('/compare_csvs', methods=['POST'])
def compare_csvs():
    """
    Compare exactly two already-uploaded CSVs (selected in the UI).
    Returns:
      {
        success: true,
        summary: { file_a: {...}, file_b: {...} },
        schema:  { shared_columns: [...], only_in_a: [...], only_in_b: [...] },
        taxon:   {
          column_used: "species_name" | "unique_name" | "otu_id" | ... | None,
          overlap_count: int, only_in_a_count: int, only_in_b_count: int,
          sample: { overlap:[..], only_in_a:[..], only_in_b:[..] }   # small examples
        }
      }
    """
    try:
        payload = request.json or {}
        names = payload.get('filenames', [])
        if not isinstance(names, list) or len(names) != 2:
            return jsonify({'success': False, 'error': 'Please provide exactly two filenames.'}), 400

        # Fetch dataframes from uploaded store
        def get_df_by_name(name):
            # exact match first
            info = chatbot.uploaded_files['csv_files'].get(name)
            if info: 
                return info.get('data'), name
            # case-insensitive fallback
            for k, v in chatbot.uploaded_files['csv_files'].items():
                if k.lower() == str(name).lower():
                    return v.get('data'), k
            return None, name

        df_a, name_a = get_df_by_name(names[0])
        df_b, name_b = get_df_by_name(names[1])

        if df_a is None:
            return jsonify({'success': False, 'error': f'Not found: {name_a}'}), 404
        if df_b is None:
            return jsonify({'success': False, 'error': f'Not found: {name_b}'}), 404

        # --- summary
        def short_label(n):
            base = os.path.basename(n)
            return os.path.splitext(base)[0]

        summary = {
            'file_a': {
                'filename': name_a,
                'short': short_label(name_a),
                'rows': int(getattr(df_a, 'shape', (0,0))[0]),
                'cols': int(getattr(df_a, 'shape', (0,0))[1]),
                'type': getattr(chatbot, 'detect_csv_type', lambda d, fn: 'generic')(df_a, name_a)
            },
            'file_b': {
                'filename': name_b,
                'short': short_label(name_b),
                'rows': int(getattr(df_b, 'shape', (0,0))[0]),
                'cols': int(getattr(df_b, 'shape', (0,0))[1]),
                'type': getattr(chatbot, 'detect_csv_type', lambda d, fn: 'generic')(df_b, name_b)
            }
        }

        # --- schema (column) comparison
        cols_a = list(df_a.columns)
        cols_b = list(df_b.columns)
        shared = sorted(list(set(cols_a) & set(cols_b)))
        only_a = sorted(list(set(cols_a) - set(cols_b)))
        only_b = sorted(list(set(cols_b) - set(cols_a)))
        schema = {
            'shared_columns': shared,
            'only_in_a': only_a,
            'only_in_b': only_b
        }

        # --- taxon/feature overlap by a best-ID column

        def pick_id_column(a_cols, b_cols):
            # preference order (kept from your original logic)
            candidates = [
                'unique_name','species_name','taxon','taxa','feature','feature_id',
                'otu_id','asv','asv_id','otu','id','name','genus','species'
            ]
            a_lower = {c.lower(): c for c in a_cols}
            b_lower = {c.lower(): c for c in b_cols}
            for cand in candidates:
                if cand in a_lower and cand in b_lower:
                    # actual names in each + normalized key
                    return a_lower[cand], b_lower[cand], cand
            return None, None, None

        def _string_set(series):

            return set(str(x).strip() for x in series.dropna().astype(str).unique())

        def auto_match_id_columns(df_a, df_b):
            """
            Generalized matcher for any pair of CSVs.
            Handles:
              - interactions (sp1/sp2) ↔ taxonomy/mapping/otu
              - otu_table (OTU_ID/ASV) ↔ mapping/taxonomy/otu_table
              - generic best string-column pair (Jaccard)
            Returns: (col_a, col_b, description) or (None, None, None)
            """
            cols_a = list(df_a.columns)
            cols_b = list(df_b.columns)
            lower_a = {c.lower(): c for c in cols_a}
            lower_b = {c.lower(): c for c in cols_b}

            # === Case 1: Interactions (sp1 & sp2) on A → match against best column in B
            if {'sp1','sp2'}.issubset(set(lower_a.keys())):
                s1 = _string_set(df_a[lower_a['sp1']])
                s2 = _string_set(df_a[lower_a['sp2']])
                idset_a = s1 | s2
                best = (0.0, None)
                for col in cols_b:
                    vals_b = _string_set(df_b[col])
                    if not vals_b:
                        continue
                    jacc = len(idset_a & vals_b) / max(1, len(idset_a | vals_b))
                    if jacc > best[0]:
                        best = (jacc, col)
                if best[1] and best[0] > 0.01:
                    return ('sp1|sp2', best[1], f'union(sp1,sp2) ↔ {best[1]}')

            # === Case 2: OTU/ASV IDs on A → match against best ID/name column in B
            for key in ['otu_id','asv','asv_id','otu','id','name']:
                if key in lower_a:
                    col_a = lower_a[key]
                    set_a = _string_set(df_a[col_a])
                    if not set_a:
                        continue
                    best = (0.0, None)
                    for col in cols_b:
                        vals_b = _string_set(df_b[col])
                        if not vals_b:
                            continue
                        jacc = len(set_a & vals_b) / max(1, len(set_a | vals_b))
                        if jacc > best[0]:
                            best = (jacc, col)
                    if best[1] and best[0] > 0.01:
                        return (col_a, best[1], f'{col_a} ↔ {best[1]}')

            # === Case 3: Generic fallback (try every string-like column pair)
            best_overlap, best_pair = 0.0, None
            for col_a in cols_a:
                vals_a = _string_set(df_a[col_a])
                if not vals_a:
                    continue
                for col_b in cols_b:
                    vals_b = _string_set(df_b[col_b])
                    if not vals_b:
                        continue
                    jacc = len(vals_a & vals_b) / max(1, len(vals_a | vals_b))
                    if jacc > best_overlap:
                        best_overlap, best_pair = jacc, (col_a, col_b)
            if best_pair and best_overlap > 0.01:
                ca, cb = best_pair
                return ca, cb, f'{ca} ↔ {cb}'

            return None, None, None

        # Try your original exact/shared-key matcher first
        col_a, col_b, _norm = pick_id_column(cols_a, cols_b)

        # If no obvious shared column, fall back to the generalized auto matcher
        if not (col_a and col_b):
            col_a, col_b, desc = auto_match_id_columns(df_a, df_b)
        else:
            desc = (col_a if col_a == col_b else f'{col_a}|{col_b}')

        # Build the taxon/feature overlap result
        taxon = {
            'column_used': None,
            'overlap_count': 0,
            'only_in_a_count': 0,
            'only_in_b_count': 0,
            'sample': {}
        }

        if col_a and col_b:
            # Support union(sp1, sp2) when requested
            if col_a == 'sp1|sp2':
                a_vals = _string_set(df_a['sp1']) | _string_set(df_a['sp2'])
            else:
                a_vals = _string_set(df_a[col_a])
            b_vals = _string_set(df_b[col_b])

            overlap = sorted(list(a_vals & b_vals))
            only_a_vals = sorted(list(a_vals - b_vals))
            only_b_vals = sorted(list(b_vals - a_vals))

            taxon.update({
                'column_used': desc,
                'overlap_count': len(overlap),
                'only_in_a_count': len(only_a_vals),
                'only_in_b_count': len(only_b_vals),
                'sample': {
                    'overlap': overlap[:10],
                    'only_in_a': only_a_vals[:10],
                    'only_in_b': only_b_vals[:10]
                }
            })
        else:
            # leave as None when no overlapping ID space exists
            taxon.update({'column_used': None})

        
        # 1. Data Quality Metrics
        data_quality = {
            'missing_values': {
                'file_a': int(df_a.isnull().sum().sum()),
                'file_b': int(df_b.isnull().sum().sum()),
                'shared_columns_missing': {col: {'file_a': int(df_a[col].isnull().sum()), 'file_b': int(df_b[col].isnull().sum())} for col in shared[:5]}
            },
            'duplicate_rows': {
                'file_a': int(df_a.duplicated().sum()),
                'file_b': int(df_b.duplicated().sum())
            },
            'data_types_match': len([c for c in shared if str(df_a[c].dtype) == str(df_b[c].dtype)]),
            'completeness_ratio': {
                'file_a': float(1 - df_a.isnull().sum().sum() / df_a.size),
                'file_b': float(1 - df_b.isnull().sum().sum() / df_b.size)
            }
        }
        
        # 2. Statistical Analysis
        numeric_shared = [
            col for col in shared
            if pd.api.types.is_numeric_dtype(df_a[col]) and pd.api.types.is_numeric_dtype(df_b[col])
        ]
        statistics = {
            'numeric_columns_count': len(numeric_shared),
            'correlations': {},
            'mean_differences': {},
            'variance_ratios': {}
        }
        
        for col in numeric_shared[:10]:  # Limit to first 10 numeric columns
            try:
                if len(df_a[col].dropna()) > 1 and len(df_b[col].dropna()) > 1:
                    # Align indices for correlation
                    common_idx = df_a.index.intersection(df_b.index)
                    if len(common_idx) > 1:
                        corr = float(df_a.loc[common_idx, col].corr(df_b.loc[common_idx, col]))
                        statistics['correlations'][col] = corr if not np.isnan(corr) else 0.0
                    
                    mean_a, mean_b = float(df_a[col].mean()), float(df_b[col].mean())
                    statistics['mean_differences'][col] = abs(mean_a - mean_b)
                    
                    var_a, var_b = float(df_a[col].var()), float(df_b[col].var())
                    if var_b != 0:
                        statistics['variance_ratios'][col] = var_a / var_b
            except:
                continue
        
       
        plots = {}

        # 1) ID overlap bar (if we found an ID column pair)
        if taxon.get('column_used'):
            plots['overlap_bar'] = {
                "data": [{
            "type": "bar",
            "x": ["Only in A", "Overlap", "Only in B"],
            "y": [
                taxon.get('only_in_a_count', 0),
                taxon.get('overlap_count', 0),
                taxon.get('only_in_b_count', 0)
            ],
            "text": [
                taxon.get('only_in_a_count', 0),
                taxon.get('overlap_count', 0),
                taxon.get('only_in_b_count', 0)
            ],
            "textposition": "auto",
            "hovertemplate": "%{x}: %{y}<extra></extra>"
        }],
        "layout": {
            "title": f"ID Overlap — {taxon.get('column_used')}",
            "yaxis": {"title": "Count"},
            "margin": {"l": 60, "r": 10, "t": 50, "b": 40}
        }
    }

        # 2) Correlation bar (from statistics)
        if statistics['correlations']:
            pairs = sorted(statistics['correlations'].items(), key=lambda x: abs(x[1]), reverse=True)
            plots['corr_bar'] = {
        "data": [{
            "type": "bar",
            "x": [p[0] for p in pairs],
            "y": [round(p[1], 4) for p in pairs],
            "hovertemplate": "%{x}: %{y}<extra></extra>"
        }],
        "layout": {
            "title": "Column Correlations (shared numeric)",
            "xaxis": {"title": "Column"},
            "yaxis": {"title": "Correlation", "range": [-1, 1]},
            "margin": {"l": 60, "r": 10, "t": 50, "b": 80}
        }
    }

        # 3) Degree histogram if File A looks like interactions (sp1/sp2)
        cols_a_lower = {c.lower() for c in df_a.columns}
        if {'sp1','sp2'}.issubset(cols_a_lower):
            s1 = df_a[[c for c in df_a.columns if c.lower() == 'sp1'][0]].astype(str)
            s2 = df_a[[c for c in df_a.columns if c.lower() == 'sp2'][0]].astype(str)
            deg = (s1.value_counts().add(s2.value_counts(), fill_value=0)).astype(int)
            plots['degree_hist'] = {
        "data": [{
            "type": "histogram",
            "x": deg.values.tolist(),
            "nbinsx": 30,
            "hovertemplate": "Degree: %{x}<extra></extra>"
        }],
        "layout": {
            "title": "Interaction Network — Node Degree Distribution (File A)",
            "xaxis": {"title": "Degree"},
            "yaxis": {"title": "Frequency"},
            "margin": {"l": 60, "r": 10, "t": 50, "b": 40}
        }
    }
        
        shared_lower = {c.lower(): c for c in shared}
        # 3. Taxonomic Hierarchy Analysis
        taxonomy_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        taxonomy_analysis = {
            'available_levels': [level for level in taxonomy_levels if level in shared],
            'level_overlap': {}
        }
        
        for level in taxonomy_levels:
            if level in shared:
                try:
                    set_a = set(str(x).strip() for x in df_a[level].dropna().astype(str).unique())
                    set_b = set(str(x).strip() for x in df_b[level].dropna().astype(str).unique())
                    overlap_vals = sorted(list(set_a & set_b))
                    unique_a_vals = sorted(list(set_a - set_b))
                    unique_b_vals = sorted(list(set_b - set_a))
                    
                    taxonomy_analysis['level_overlap'][level] = {
                        'shared_count': len(overlap_vals),
                        'unique_a_count': len(unique_a_vals),
                        'unique_b_count': len(unique_b_vals),
                        'jaccard_similarity': len(overlap_vals) / len(set_a | set_b) if len(set_a | set_b) > 0 else 0,
                        'sample_shared': overlap_vals[:5],
                        'sample_unique_a': unique_a_vals[:5],
                        'sample_unique_b': unique_b_vals[:5]
                    }
                except:
                    continue
        
        # 4. Distribution Analysis
        numeric_cols_a = df_a.select_dtypes(include=[np.number]).columns
        numeric_cols_b = df_b.select_dtypes(include=[np.number]).columns
        
        distribution_analysis = {
            'abundance_ranges': {
                'file_a': {
                    'min': float(df_a[numeric_cols_a].min().min()) if len(numeric_cols_a) > 0 else 0,
                    'max': float(df_a[numeric_cols_a].max().max()) if len(numeric_cols_a) > 0 else 0,
                    'median': float(df_a[numeric_cols_a].median().median()) if len(numeric_cols_a) > 0 else 0
                },
                'file_b': {
                    'min': float(df_b[numeric_cols_b].min().min()) if len(numeric_cols_b) > 0 else 0,
                    'max': float(df_b[numeric_cols_b].max().max()) if len(numeric_cols_b) > 0 else 0,
                    'median': float(df_b[numeric_cols_b].median().median()) if len(numeric_cols_b) > 0 else 0
                }
            },
            'zero_abundance_ratio': {
                'file_a': float((df_a[numeric_cols_a] == 0).sum().sum() / df_a[numeric_cols_a].size) if len(numeric_cols_a) > 0 and df_a[numeric_cols_a].size > 0 else 0,
                'file_b': float((df_b[numeric_cols_b] == 0).sum().sum() / df_b[numeric_cols_b].size) if len(numeric_cols_b) > 0 and df_b[numeric_cols_b].size > 0 else 0
            },
            'sparsity_comparison': {
                'file_a_sparse': float((df_a[numeric_cols_a] == 0).sum().sum() / df_a[numeric_cols_a].size > 0.5) if len(numeric_cols_a) > 0 else False,
                'file_b_sparse': float((df_b[numeric_cols_b] == 0).sum().sum() / df_b[numeric_cols_b].size > 0.5) if len(numeric_cols_b) > 0 else False
            }
        }
        
        # 5. Sample/Column Similarity
        jaccard_sim = len(set(df_a.columns) & set(df_b.columns)) / len(set(df_a.columns) | set(df_b.columns)) if len(set(df_a.columns) | set(df_b.columns)) > 0 else 0
        
        sample_analysis = {
            'column_overlap_count': len(shared),
            'jaccard_similarity': float(jaccard_sim),
            'similarity_category': 'high' if jaccard_sim > 0.7 else 'medium' if jaccard_sim > 0.4 else 'low',
            'total_unique_columns': len(set(df_a.columns) | set(df_b.columns)),
            'column_overlap_percentage': float(len(shared) / len(set(df_a.columns) | set(df_b.columns)) * 100) if len(set(df_a.columns) | set(df_b.columns)) > 0 else 0
        }
        
        # 6. Diversity Metrics (for microbiome data)
        def calculate_shannon_diversity(df_numeric):
            if df_numeric.empty:
                return 0
            # Calculate Shannon diversity for each sample (column)
            diversities = []
            for col in df_numeric.columns:
                abundances = df_numeric[col].values
                abundances = abundances[abundances > 0]  # Remove zeros
                if len(abundances) > 0:
                    proportions = abundances / abundances.sum()
                    shannon = -np.sum(proportions * np.log(proportions))
                    diversities.append(shannon)
            return float(np.mean(diversities)) if diversities else 0
        
        diversity_metrics = {
            'shannon_diversity': {
                'file_a': calculate_shannon_diversity(df_a[numeric_cols_a]) if len(numeric_cols_a) > 0 else 0,
                'file_b': calculate_shannon_diversity(df_b[numeric_cols_b]) if len(numeric_cols_b) > 0 else 0
            },
            'richness': {
                'file_a': int((df_a[numeric_cols_a] > 0).sum().sum()) if len(numeric_cols_a) > 0 else 0,
                'file_b': int((df_b[numeric_cols_b] > 0).sum().sum()) if len(numeric_cols_b) > 0 else 0
            },
            'total_features': {
                'file_a': len(df_a),
                'file_b': len(df_b)
            }
        }
        
        # 7. Visualization Data
        visualization_data = {
            'venn_diagram': {
                'overlap': len(overlap) if 'overlap' in locals() else 0,
                'only_a': len(only_a_vals) if 'only_a_vals' in locals() else 0,
                'only_b': len(only_b_vals) if 'only_b_vals' in locals() else 0,
                'total_unique': len((a_vals | b_vals)) if 'a_vals' in locals() and 'b_vals' in locals() else 0
            },
            'correlation_summary': {
                'high_correlation_count': len([v for v in statistics['correlations'].values() if abs(v) > 0.7]),
                'medium_correlation_count': len([v for v in statistics['correlations'].values() if 0.3 < abs(v) <= 0.7]),
                'low_correlation_count': len([v for v in statistics['correlations'].values() if abs(v) <= 0.3])
            }
        }
        
        # 8. Smart Recommendations
        recommendations = {
            'data_compatibility': sample_analysis['similarity_category'],
            'merge_feasibility': 'feasible' if jaccard_sim > 0.5 else 'challenging' if jaccard_sim > 0.2 else 'difficult',
            'suggested_actions': [],
            'quality_score': float((data_quality['completeness_ratio']['file_a'] + data_quality['completeness_ratio']['file_b']) / 2)
        }
        
        # Generate specific recommendations
        if data_quality['missing_values']['file_a'] > df_a.size * 0.1 or data_quality['missing_values']['file_b'] > df_b.size * 0.1:
            recommendations['suggested_actions'].append('Handle missing values before analysis')
        
        # Check for division by zero before comparing abundance ranges
        max_a = distribution_analysis['abundance_ranges']['file_a']['max']
        max_b = distribution_analysis['abundance_ranges']['file_b']['max']
        if max_a > 0 and max_b > 0:
            if max_a / max_b > 10 or max_b / max_a > 10:
                recommendations['suggested_actions'].append('Consider normalizing abundance values')
        
        if len(statistics['correlations']) > 0 and np.mean(list(statistics['correlations'].values())) < 0.3:
            recommendations['suggested_actions'].append('Check for batch effects or systematic differences')
        
        if jaccard_sim < 0.3:
            recommendations['suggested_actions'].append('Files have low overlap - verify they are comparable datasets')
        
        if not recommendations['suggested_actions']:
            recommendations['suggested_actions'].append('Files appear compatible for analysis')


        def _numeric_cols(df):
            return [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]

        def _as_comp_sorted(values: pd.Series) -> np.ndarray:
            """Coerce to numeric, clip negatives to 0, normalize to composition, sort desc."""
            arr = pd.to_numeric(values, errors='coerce').fillna(0.0).to_numpy(dtype=float)
            # no negatives in counts
            arr = np.clip(arr, 0, None)
            s = arr.sum()
            if s > 0:
                arr = arr / s
            # sorts descending for an ID-agnostic rank-abundance comparison
            return np.sort(arr)[::-1]

        def _pad_to_same_len(a: np.ndarray, b: np.ndarray):
            if len(a) == len(b): return a, b
            if len(a) < len(b):
                a = np.pad(a, (0, len(b)-len(a)), mode='constant')
            else:
                b = np.pad(b, (0, len(a)-len(b)), mode='constant')
            return a, b

        def _spearman_from_arrays(x: np.ndarray, y: np.ndarray):
            if len(x) < 2 or len(y) < 2: return float('nan')
            # Spearman via rank + Pearson
            rx = pd.Series(x).rank(method='average')
            ry = pd.Series(y).rank(method='average')
            rho = rx.corr(ry)
            return float(rho) if pd.notna(rho) else float('nan')

        numeric_fallback = None
        # Only invoke fallback if we **did not** find a usable ID column pair
        if not taxon.get('column_used'):
            num_a = _numeric_cols(df_a)
            num_b = _numeric_cols(df_b)
            shared_samples = [c for c in num_a if c in num_b]

            per_sample = []
            for col in shared_samples:
                xa = _as_comp_sorted(df_a[col])
                xb = _as_comp_sorted(df_b[col])
                xa, xb = _pad_to_same_len(xa, xb)

                # Bray–Curtis dissimilarity 
                # Similarity = 1 - D (It is nicer for bar plots where higher is better)
                bc_diss = float(0.5 * np.abs(xa - xb).sum())
                bc_sim  = float(1.0 - bc_diss)

                rho = _spearman_from_arrays(xa, xb)
                # Pearson on shapes 
                pear = float(pd.Series(xa).corr(pd.Series(xb))) if len(xa) > 1 else float('nan')

                per_sample.append({
                    "sample": col,
                    "bray_curtis_dissimilarity": bc_diss,
                    "bray_curtis_similarity": bc_sim,
                    "spearman": rho,
                    "pearson": pear
                })

            if per_sample:
                numeric_fallback = {
                    "shared_samples": shared_samples,
                    "per_sample": per_sample,
                    "summary": {
                        "mean_bray_curtis_dissimilarity": float(np.mean([d["bray_curtis_dissimilarity"] for d in per_sample])),
                        "mean_bray_curtis_similarity":   float(np.mean([d["bray_curtis_similarity"]   for d in per_sample])),
                        "mean_spearman":                 float(np.nanmean([d["spearman"]               for d in per_sample])),
                        "mean_pearson":                  float(np.nanmean([d["pearson"]                for d in per_sample]))
                    }
                }

                
                samples   = [d["sample"] for d in per_sample]
                bc_sim    = [round(d["bray_curtis_similarity"], 6) for d in per_sample]
                rho_vals  = [None if np.isnan(d["spearman"]) else round(d["spearman"], 6) for d in per_sample]

                plots['per_sample_bray'] = {
                    "data": [{
                        "type": "bar",
                        "x": samples,
                        "y": bc_sim,
                        "hovertemplate": "%{x}: %{y:.3f}<extra></extra>"
                    }],
                    "layout": {
                        "title": "Per-sample Similarity (1 − Bray–Curtis) — no row ID match",
                        "xaxis": {"title": "Sample", "automargin": True},
                        "yaxis": {"title": "Similarity (higher = more similar)", "range": [0, 1]},
                        "margin": {"l": 60, "r": 10, "t": 50, "b": 80}
                    }
                }

                plots['per_sample_spearman'] = {
                    "data": [{
                        "type": "bar",
                        "x": samples,
                        "y": rho_vals,
                        "hovertemplate": "%{x}: ρ=%{y:.3f}<extra></extra>"
                    }],
                    "layout": {
                        "title": "Per-sample Spearman (rank-abundance) — no row ID match",
                        "xaxis": {"title": "Sample", "automargin": True},
                        "yaxis": {"title": "ρ", "range": [-1, 1]},
                        "margin": {"l": 60, "r": 10, "t": 50, "b": 80}
                    }
                }

        return jsonify({
            'success': True, 
            'summary': summary, 
            'schema': schema, 
            'taxon': taxon,
            'data_quality': data_quality,
            'statistics': statistics,
            'taxonomy_analysis': taxonomy_analysis,
            'distribution_analysis': distribution_analysis,
            'sample_analysis': sample_analysis,
            'diversity_metrics': diversity_metrics,
            'visualization_data': visualization_data,
            'recommendations': recommendations,
            'numeric_fallback': numeric_fallback,
            'plots': plots
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 400


@app.route('/history')
def history():
    """Get chat history"""
    return jsonify({'success': True, 'chat_history': chatbot.chat_history})

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

        # Registers the uploaded file in-memory so depth queries work immediately
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
            result = chatbot.analyze_csv_file(file_path)
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
                    continue  # Skips non-PDF files
                
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


@app.route('/health', methods=['GET'])
def health():
    exists = {
        'DATA_DIR': os.path.isdir(getattr(chatbot, 'data_dir', 'data')),
        'UPLOAD_DIR': os.path.isdir(getattr(chatbot, 'upload_dir', 'uploads'))
    }
    return jsonify({'ok': True, 'exists': exists})

@app.route('/docs', methods=['GET'])
def docs():
    return jsonify({
        "app": "Microbiome XAI Chatbot (Pro v3)",
        "endpoints": [
            "/", "/chat", "/upload_file", "/get_uploaded_files",
            "/history", "/clear_history",
            "/compare_interactions", "/compare_taxonomy", "/compare_mappings",
            "/health", "/docs"
        ]
    })

@app.route('/compare_interactions', methods=['GET'])
def compare_interactions_route():
    a = request.args.get('a', 'interactions_enhanced')
    b = request.args.get('b', 'interactions_short')
    a_df, a_name = chatbot._load_csv_from_known_locations(a)
    b_df, b_name = chatbot._load_csv_from_known_locations(b)
    if a_df is None or b_df is None:
        return jsonify({"success": False, "error": f"Could not find both CSVs: a='{a}', b='{b}'"}), 400
    res = chatbot.compare_interactions(a_df, b_df)
    return jsonify({"success": True, "a": a_name, "b": b_name, **res})

@app.route('/compare_taxonomy', methods=['GET'])
def compare_taxonomy_route():
    a = request.args.get('a', 'taxonomy')
    b = request.args.get('b', 'taxonomy_table')
    a_df, a_name = chatbot._load_csv_from_known_locations(a)
    b_df, b_name = chatbot._load_csv_from_known_locations(b)
    if a_df is None or b_df is None:
        return jsonify({"success": False, "error": f"Could not find two taxonomy CSVs to compare (a='{a}', b='{b}')."}), 400
    res = chatbot.compare_taxonomy(a_df, b_df)
    return jsonify({"success": True, "a": a_name, "b": b_name, **res})

@app.route('/compare_mappings', methods=['GET'])
def compare_mappings_route():
    a = request.args.get('a', 'mapping')
    b = request.args.get('b', 'augmented')
    a_df, a_name = chatbot._load_csv_from_known_locations(a)
    b_df, b_name = chatbot._load_csv_from_known_locations(b)
    if a_df is None or b_df is None:
        return jsonify({"success": False, "error": f"Could not find two mapping CSVs to compare (a='{a}', b='{b}')."}), 400
    res = chatbot.compare_mappings(a_df, b_df)
    return jsonify({"success": True, "a": a_name, "b": b_name, **res})

@app.route('/explain_detailed', methods=['POST'])
def explain_detailed():
    data = request.get_json(force=True)
    s1 = data.get('species1'); s2 = data.get('species2'); inferred = data.get('label', 'Unknown')
    trace = chatbot.proof_trace_for_interaction(s1, s2, inferred)
    return jsonify({'success': 'error' not in trace, **trace})

@app.route('/network_null_test', methods=['GET'])
def network_null_test():
    out = chatbot.network_null_significance(num_draws=int(request.args.get('n', 200)))
    return jsonify({'success': 'error' not in out, **out})

@app.route('/delete_csvs', methods=['POST'])
def delete_csvs():
    """Delete selected CSV files"""
    try:
        data = request.json
        filenames = data.get('filenames', [])
        
        if not filenames:
            return jsonify({'success': False, 'error': 'No filenames provided'}), 400
        
        deleted_files = []
        errors = []
        
        for filename in filenames:
            try:
                # Remove from uploaded_files
                if filename in chatbot.uploaded_files['csv_files']:
                    del chatbot.uploaded_files['csv_files'][filename]
                    deleted_files.append(filename)
                else:
                    errors.append(f"File '{filename}' not found")
            except Exception as e:
                errors.append(f"Error deleting '{filename}': {str(e)}")
        
        return jsonify({
            'success': True,
            'deleted_files': deleted_files,
            'errors': errors,
            'message': f"Deleted {len(deleted_files)} file(s)"
        })
        
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.get("/chat_help")
def chat_help():
    return jsonify({"ok": True, "message": "Chatbot online."})

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=5000, threaded=False, use_reloader=False)