#!/usr/bin/env python3
"""
RNA-seq Differential Expression Analysis Pipeline
This script processes GEO2R output, performs pathway analysis using DAVID,
and identifies drugs from KEGG that target disease-related genes.
"""

import pandas as pd
import requests
import time
import argparse
import sys
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION SECTION - YOU MAY NEED TO MODIFY THESE PARAMETERS
# ============================================================================

class Config:
    """Configuration class for pipeline parameters"""
    
    # File paths - MODIFY THESE AS NEEDED
    INPUT_FILE = "input.tsv"  # Change to your input file path
    OUTPUT_FILE = "rna_seq_analysis_results.xlsx"  # Change to desired output file name
    
    # Significance thresholds - MODIFY IF NEEDED
    P_VALUE_THRESHOLD = 0.05  # Standard p-value threshold
    LOGFC_THRESHOLD = 2.0  # |logFC| > 2 as specified
    
    # Column names - These will be detected dynamically, but you can specify if known
    # If your TSV has different column names, the script will try to detect them
    # But you can manually specify them here if detection fails
    GENE_COLUMN_POSSIBLE_NAMES = ['Gene.symbol', 'Gene', 'Symbol', 'gene', 'GeneSymbol', 'Gene Name']
    P_VALUE_COLUMN_POSSIBLE_NAMES = ['P.Value', 'pvalue', 'PValue', 'p-value', 'p_val', 'P', 'p']
    LOGFC_COLUMN_POSSIBLE_NAMES = ['logFC', 'LogFC', 'log2FC', 'log2FoldChange', 'FoldChange', 'log_fc']
    
    # Enrichr API endpoints (if you want to use Enrichr instead of DAVID, you can modify the code to use these)
    ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
    ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"
    
    KEGG_BASE_URL = "http://rest.kegg.jp"
    
    # Rate limiting (be respectful to APIs)
    API_DELAY = 0.5  # Delay between API calls in seconds

# ============================================================================
# MAIN PIPELINE CLASS
# ============================================================================

class RNASeqAnalysisPipeline:
    """Main pipeline class for RNA-seq analysis"""
    
    def __init__(self, disease_name: str, input_file: str = None):
        """
        Initialize the pipeline with disease name and optional input file
        
        Args:
            disease_name: Name of the disease to analyze
            input_file: Path to input TSV file (overrides Config.INPUT_FILE if provided)
        """
        self.disease_name = disease_name
        self.input_file = input_file if input_file else Config.INPUT_FILE
        self.df = None
        self.upregulated_genes = []
        self.downregulated_genes = []
        self.pathway_results = {}
        self.drug_results = {}
        
        # Validate API credentials
        '''if Config.DAVID_EMAIL == "your.email@example.com" or Config.DAVID_API_KEY == "YOUR_DAVID_API_KEY_HERE":
            print("WARNING: Please set your DAVID email and API key in the Config class")
            print("You can register for an API key at: https://david.ncifcrf.gov/content.jsp?file=WS.html")'''
        
    def load_and_validate_input(self) -> pd.DataFrame:
        """
        Load the TSV file and dynamically detect column names
        
        Returns:
            DataFrame with the loaded data
        """
        print(f"\n{'='*60}")
        print(f"Step 1: Loading input file: {self.input_file}")
        print(f"{'='*60}")
        
        try:
            # Load the TSV file
            self.df = pd.read_csv(self.input_file, sep='\t')
            print(f"File loaded successfully. Shape: {self.df.shape}")
            print(f"Columns found: {list(self.df.columns)}")
            
            # Dynamically detect column names
            self._detect_columns()
            
            # Filter out rows with missing gene symbols
            initial_rows = len(self.df)
            self.df = self.df.dropna(subset=[self.gene_column])
            self.df = self.df[self.df[self.gene_column].astype(str) != 'nan']
            self.df = self.df[self.df[self.gene_column].astype(str) != '']
            print(f"Removed {initial_rows - len(self.df)} rows with missing gene symbols")
            
            return self.df
            
        except FileNotFoundError:
            print(f"ERROR: Input file '{self.input_file}' not found!")
            print("Please check the file path or set Config.INPUT_FILE correctly")
            sys.exit(1)
        except Exception as e:
            print(f"ERROR loading file: {e}")
            sys.exit(1)
    
    def _detect_columns(self):
        """
        Dynamically detect gene symbol, p-value, and logFC columns
        """
        # Detect gene column
        self.gene_column = None
        for col in self.df.columns:
            if any(gene_name in col for gene_name in Config.GENE_COLUMN_POSSIBLE_NAMES):
                self.gene_column = col
                break
        
        if not self.gene_column:
            # If no match found, use first string column
            for col in self.df.columns:
                if self.df[col].dtype == 'object':
                    self.gene_column = col
                    print(f"Using '{col}' as gene column (auto-detected)")
                    break
        
        if not self.gene_column:
            raise ValueError("Could not detect gene symbol column")
        
        # Detect p-value column
        self.pval_column = None
        for col in self.df.columns:
            if any(pval_name in col for pval_name in Config.P_VALUE_COLUMN_POSSIBLE_NAMES):
                self.pval_column = col
                break
        
        if not self.pval_column:
            # Try to find numeric column with values between 0 and 1
            for col in self.df.columns:
                if self.df[col].dtype in ['float64', 'int64']:
                    if ((self.df[col] >= 0) & (self.df[col] <= 1)).mean() > 0.5:
                        self.pval_column = col
                        break
        
        if not self.pval_column:
            raise ValueError("Could not detect p-value column")
        
        # Detect logFC column
        self.logfc_column = None
        for col in self.df.columns:
            if any(logfc_name in col for logfc_name in Config.LOGFC_COLUMN_POSSIBLE_NAMES):
                self.logfc_column = col
                break
        
        if not self.logfc_column:
            # Try to find numeric column with both positive and negative values
            for col in self.df.columns:
                if self.df[col].dtype in ['float64', 'int64'] and col != self.pval_column:
                    if (self.df[col] > 0).any() and (self.df[col] < 0).any():
                        self.logfc_column = col
                        break
        
        if not self.logfc_column:
            raise ValueError("Could not detect logFC column")
        
        print(f"\nDetected columns:")
        print(f"  Gene column: {self.gene_column}")
        print(f"  P-value column: {self.pval_column}")
        print(f"  LogFC column: {self.logfc_column}")
    
    def filter_significant_genes(self) -> Tuple[List[str], List[str]]:
        """
        Filter genes based on p-value and logFC thresholds
        
        Returns:
            Tuple of (upregulated_genes, downregulated_genes)
        """
        print(f"\n{'='*60}")
        print(f"Step 2: Filtering significant genes")
        print(f"Thresholds: p-value < {Config.P_VALUE_THRESHOLD}, |logFC| > {Config.LOGFC_THRESHOLD}")
        print(f"{'='*60}")
        
        # Create significance mask
        sig_mask = (self.df[self.pval_column] < Config.P_VALUE_THRESHOLD) & \
                   (abs(self.df[self.logfc_column]) > Config.LOGFC_THRESHOLD)
        
        significant_df = self.df[sig_mask].copy()
        print(f"Found {len(significant_df)} significant genes")
        
        # Separate up and down regulated
        up_mask = significant_df[self.logfc_column] > Config.LOGFC_THRESHOLD
        self.upregulated_genes = significant_df[up_mask][self.gene_column].tolist()
        self.downregulated_genes = significant_df[~up_mask][self.gene_column].tolist()
        
        # Create DataFrames for up and down regulated
        self.up_df = significant_df[up_mask][[self.gene_column, self.pval_column, self.logfc_column]].copy()
        self.up_df['Regulation'] = 'Upregulated'
        
        self.down_df = significant_df[~up_mask][[self.gene_column, self.pval_column, self.logfc_column]].copy()
        self.down_df['Regulation'] = 'Downregulated'
        
        print(f"Upregulated genes: {len(self.upregulated_genes)}")
        print(f"Downregulated genes: {len(self.downregulated_genes)}")
        
        return self.upregulated_genes, self.downregulated_genes
    
    def query_enrichr_pathways(self, gene_list: List[str]) -> Dict:

        print("\n============================================================")
        print("Step 3: Querying Enrichr for pathway analysis")
        print("============================================================")

        genes = "\n".join(gene_list)

        payload = {
            'list': (None, genes),
            'description': (None, 'RNAseq analysis')
        }

        response = requests.post(Config.ENRICHR_ADD_URL, files=payload)
        if not response.ok:
            print("Error uploading gene list to Enrichr")
            return {}

        user_list_id = response.json()['userListId']

        params = {
            'userListId': user_list_id,
            'backgroundType': 'KEGG_2021_Human'
        }

        response = requests.get(Config.ENRICHR_ENRICH_URL, params=params)

        pathways = {}

        if response.ok:
            data = response.json()['KEGG_2021_Human']

            for row in data[:20]:
                pathway_name = row[1]
                p_value = row[2]
                genes = row[5]

                pathways[pathway_name] = {
                    "p_value": p_value,
                    "genes": genes[:10]
                }

        print(f"Found {len(pathways)} enriched pathways")

        return pathways
    
    def search_kegg_disease(self) -> Dict:
        """
        Search KEGG DISEASE database for disease-related genes
        
        Returns:
            Dictionary with disease-related genes and information
        """
        print(f"\n{'='*60}")
        print(f"Step 4: Searching KEGG DISEASE for: {self.disease_name}")
        print(f"{'='*60}")
        
        disease_genes = {}
        
        try:
            # Search for disease in KEGG
            search_url = f"{Config.KEGG_BASE_URL}/find/disease/{self.disease_name}"
            response = requests.get(search_url)
            time.sleep(Config.API_DELAY)
            
            if response.status_code == 200 and response.text.strip():
                lines = response.text.strip().split('\n')
                disease_ids = []
                
                # Parse disease IDs from search results
                for line in lines[:5]:  # Take top 5 matches
                    parts = line.split('\t')
                    if parts:
                        disease_id = parts[0]
                        disease_ids.append(disease_id)
                        print(f"Found disease entry: {line[:100]}...")
                
                # Get gene associations for each disease ID
                for disease_id in disease_ids:
                    disease_url = f"{Config.KEGG_BASE_URL}/get/{disease_id}"
                    response = requests.get(disease_url)
                    time.sleep(Config.API_DELAY)
                    
                    if response.status_code == 200:
                        content = response.text
                        # Extract gene information
                        genes = []
                        in_gene_section = False
                        
                        for line in content.split('\n'):
                            if line.startswith('GENE'):
                                in_gene_section = True
                                # Parse gene line
                                parts = line.split()
                                if len(parts) >= 2:
                                    gene_symbol = parts[1].split('(')[0].strip()
                                    genes.append(gene_symbol)
                            elif in_gene_section and line.startswith(' '):
                                # Continuation of gene section
                                parts = line.strip().split()
                                if parts:
                                    gene_symbol = parts[0].split('(')[0].strip()
                                    genes.append(gene_symbol)
                            elif line.strip() and not line.startswith(' ') and in_gene_section:
                                in_gene_section = False
                        
                        disease_genes[disease_id] = {
                            'name': disease_id.split(':')[-1] if ':' in disease_id else disease_id,
                            'genes': genes[:50]  # Limit to first 50 genes
                        }
                        
                        print(f"  {disease_id}: Found {len(genes)} associated genes")
            
            return disease_genes
            
        except Exception as e:
            print(f"Error searching KEGG DISEASE: {e}")
            return {}
    
    def get_drug_targets(self, gene_list: List[str]) -> Dict:
        """
        Query KEGG DRUG for drugs targeting specific genes
        
        Args:
            gene_list: List of gene symbols
        
        Returns:
            Dictionary with drug information for each gene
        """
        print(f"\n{'='*60}")
        print(f"Step 5: Finding drugs targeting disease genes")
        print(f"{'='*60}")
        
        drug_info = {}
        
        for i, gene in enumerate(gene_list[:50]):  # Limit to first 50 genes
            if i % 10 == 0:
                print(f"Processing gene {i+1}/{min(50, len(gene_list))}: {gene}")
            
            try:
                # Search for drugs targeting this gene
                # KEGG DRUG search by target
                search_url = f"{Config.KEGG_BASE_URL}/find/drug/{gene}"
                response = requests.get(search_url)
                time.sleep(Config.API_DELAY)
                
                if response.status_code == 200 and response.text.strip():
                    drugs = []
                    lines = response.text.strip().split('\n')
                    
                    for line in lines[:5]:  # Limit to top 5 drugs per gene
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            drug_id = parts[0]
                            drug_name = parts[1]
                            
                            # Get detailed drug info
                            drug_url = f"{Config.KEGG_BASE_URL}/get/{drug_id}"
                            drug_response = requests.get(drug_url)
                            time.sleep(Config.API_DELAY)
                            
                            if drug_response.status_code == 200:
                                drug_content = drug_response.text
                                drug_info_text = self._parse_drug_info(drug_content, gene)
                                if drug_info_text:
                                    drugs.append({
                                        'drug_id': drug_id,
                                        'drug_name': drug_name,
                                        'info': drug_info_text
                                    })
                    
                    if drugs:
                        drug_info[gene] = drugs
                        print(f"  Found {len(drugs)} drugs for {gene}")
                        
            except Exception as e:
                print(f"  Error processing gene {gene}: {e}")
                continue
        
        return drug_info
    
    def _parse_drug_info(self, drug_content: str, target_gene: str) -> Optional[str]:
        """
        Parse KEGG drug information to extract relevant details
        
        Args:
            drug_content: Raw KEGG drug entry content
            target_gene: Target gene name
        
        Returns:
            Formatted drug information string
        """
        info_lines = []
        in_target_section = False
        in_comment_section = False
        
        for line in drug_content.split('\n'):
            if line.startswith('TARGET'):
                in_target_section = True
                if target_gene.lower() in line.lower():
                    info_lines.append(line.strip())
            elif in_target_section and line.startswith(' '):
                if target_gene.lower() in line.lower():
                    info_lines.append(line.strip())
            elif line.strip() and not line.startswith(' ') and in_target_section:
                in_target_section = False
            
            if line.startswith('COMMENT'):
                in_comment_section = True
                info_lines.append(line.strip())
            elif in_comment_section and line.startswith(' '):
                info_lines.append(line.strip())
            elif line.strip() and not line.startswith(' ') and in_comment_section:
                in_comment_section = False
        
        return '\n'.join(info_lines) if info_lines else None
    
    def create_output_excel(self):
        """
        Create Excel output with multiple sheets
        """
        print(f"\n{'='*60}")
        print(f"Step 6: Creating Excel output: {Config.OUTPUT_FILE}")
        print(f"{'='*60}")
        
        with pd.ExcelWriter(Config.OUTPUT_FILE, engine='openpyxl') as writer:
            # Sheet 1: Original input data
            self.df.to_excel(writer, sheet_name='Original_Data', index=False)
            
            # Sheet 2: Upregulated genes
            if hasattr(self, 'up_df') and not self.up_df.empty:
                self.up_df.to_excel(writer, sheet_name='Upregulated_Genes', index=False)
            else:
                pd.DataFrame({'Message': ['No upregulated genes found']}).to_excel(
                    writer, sheet_name='Upregulated_Genes', index=False)
            
            # Sheet 3: Downregulated genes
            if hasattr(self, 'down_df') and not self.down_df.empty:
                self.down_df.to_excel(writer, sheet_name='Downregulated_Genes', index=False)
            else:
                pd.DataFrame({'Message': ['No downregulated genes found']}).to_excel(
                    writer, sheet_name='Downregulated_Genes', index=False)
            
            # Sheet 4: Pathway analysis results
            pathway_data = []
            for pathway, info in self.pathway_results.items():
                pathway_data.append({
                    'Pathway': pathway,
                    'P_Value': info.get('p_value', 'N/A'),
                    'Genes': ', '.join(info.get('genes', [])[:10])
                })
            
            if pathway_data:
                pd.DataFrame(pathway_data).to_excel(writer, sheet_name='Pathways', index=False)
            else:
                pd.DataFrame({'Message': ['No pathway data available']}).to_excel(
                    writer, sheet_name='Pathways', index=False)
            
            # Sheet 5: Drug information
            drug_data = []
            for gene, drugs in self.drug_results.items():
                for drug in drugs[:3]:  # Limit to top 3 drugs per gene
                    drug_data.append({
                        'Gene': gene,
                        'Drug_ID': drug.get('drug_id', ''),
                        'Drug_Name': drug.get('drug_name', ''),
                        'Drug_Info': drug.get('info', '')[:500]  # Truncate long text
                    })
            
            if drug_data:
                pd.DataFrame(drug_data).to_excel(writer, sheet_name='Drugs', index=False)
            else:
                pd.DataFrame({'Message': ['No drug data available']}).to_excel(
                    writer, sheet_name='Drugs', index=False)
            
            # Sheet 6: Final combined output
            final_data = []
            
            # Combine all significant genes
            all_sig_genes = []
            if hasattr(self, 'up_df') and not self.up_df.empty:
                all_sig_genes.extend(self.up_df.to_dict('records'))
            if hasattr(self, 'down_df') and not self.down_df.empty:
                all_sig_genes.extend(self.down_df.to_dict('records'))
            
            for gene_info in all_sig_genes:
                gene_symbol = gene_info[self.gene_column]
                
                # Find pathways for this gene
                gene_pathways = []
                for pathway, info in self.pathway_results.items():
                    if gene_symbol in info.get('genes', []):
                        gene_pathways.append(pathway)
                
                # Find drugs for this gene
                gene_drugs = []
                if gene_symbol in self.drug_results:
                    for drug in self.drug_results[gene_symbol][:3]:
                        gene_drugs.append(f"{drug.get('drug_name', '')} ({drug.get('drug_id', '')})")
                
                final_data.append({
                    'Gene': gene_symbol,
                    'P_Value': gene_info[self.pval_column],
                    'LogFC': gene_info[self.logfc_column],
                    'Regulation': gene_info.get('Regulation', 'Unknown'),
                    'Pathways': '; '.join(gene_pathways[:5]) if gene_pathways else 'No pathways found',
                    'Drugs': '; '.join(gene_drugs) if gene_drugs else 'No drugs found',
                    'Disease': self.disease_name
                })
            
            if final_data:
                pd.DataFrame(final_data).to_excel(writer, sheet_name='Final_Results', index=False)
            else:
                pd.DataFrame({'Message': ['No significant genes found']}).to_excel(
                    writer, sheet_name='Final_Results', index=False)
        
        print(f"Excel file created successfully: {Config.OUTPUT_FILE}")
    
    def run_pipeline(self):
        """
        Execute the complete pipeline
        """
        print(f"\n{'='*60}")
        print(f"STARTING RNA-SEQ ANALYSIS PIPELINE")
        print(f"Disease: {self.disease_name}")
        print(f"{'='*60}")
        
        # Step 1: Load data
        self.load_and_validate_input()
        
        # Step 2: Filter significant genes
        up_genes, down_genes = self.filter_significant_genes()
        
        # Step 3: Pathway analysis (combine up and down regulated)
        all_sig_genes = up_genes + down_genes
        if all_sig_genes:
            self.pathway_results = self.query_enrichr_pathways(all_sig_genes)
        
        # Step 4: Search KEGG disease
        disease_genes_info = self.search_kegg_disease()
        
        # Find overlap between significant genes and disease genes
        disease_related_genes = []
        for disease_id, info in disease_genes_info.items():
            overlap = set(all_sig_genes) & set(info.get('genes', []))
            if overlap:
                disease_related_genes.extend(list(overlap))
                print(f"\nOverlap with {disease_id}: {len(overlap)} genes")
                print(f"  Example genes: {list(overlap)[:5]}")
        
        disease_related_genes = list(set(disease_related_genes))
        
        # Step 5: Find drugs for disease-related genes
        if disease_related_genes:
            self.drug_results = self.get_drug_targets(disease_related_genes)
        
        # Step 6: Create output
        self.create_output_excel()
        
        print(f"\n{'='*60}")
        print(f"PIPELINE COMPLETED SUCCESSFULLY")
        print(f"{'='*60}")
        print(f"Summary:")
        print(f"  Total genes analyzed: {len(self.df)}")
        print(f"  Significant genes: {len(up_genes) + len(down_genes)}")
        print(f"    Upregulated: {len(up_genes)}")
        print(f"    Downregulated: {len(down_genes)}")
        print(f"  Pathways found: {len(self.pathway_results)}")
        print(f"  Disease-related genes: {len(disease_related_genes)}")
        print(f"  Genes with drug targets: {len(self.drug_results)}")
        print(f"\nResults saved to: {Config.OUTPUT_FILE}")

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def main():
    """Main function to run the pipeline from command line"""
    parser = argparse.ArgumentParser(
        description='RNA-seq Analysis Pipeline for Disease Gene Discovery',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python rna_seq_pipeline.py "Alzheimer disease"
  python rna_seq_pipeline.py "Parkinson disease" --input my_data.tsv --output results.xlsx
  python rna_seq_pipeline.py "breast cancer" --pval 0.01 --logfc 1.5
        """
    )
    
    parser.add_argument(
        'disease',
        type=str,
        help='Name of the disease to analyze (e.g., "Alzheimer disease", "colorectal cancer ")'
    )
    
    parser.add_argument(
        '--input', '-i',
        type=str,
        default=Config.INPUT_FILE,
        help=f'Input TSV file path (default: {Config.INPUT_FILE})'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=Config.OUTPUT_FILE,
        help=f'Output Excel file path (default: {Config.OUTPUT_FILE})'
    )
    
    parser.add_argument(
        '--pval',
        type=float,
        default=Config.P_VALUE_THRESHOLD,
        help=f'P-value threshold (default: {Config.P_VALUE_THRESHOLD})'
    )
    
    parser.add_argument(
        '--logfc',
        type=float,
        default=Config.LOGFC_THRESHOLD,
        help=f'LogFC threshold (default: {Config.LOGFC_THRESHOLD})'
    )
    
    '''parser.add_argument(
        '--email',
        type=str,
        default=Config.DAVID_EMAIL,
        help='Your email for DAVID API registration'
    )
    
    parser.add_argument(
        '--apikey',
        type=str,
        default=Config.DAVID_API_KEY,
        help='Your DAVID API key'
    )'''
    
    args = parser.parse_args()
    
    # Update configuration with command line arguments
    Config.INPUT_FILE = args.input
    Config.OUTPUT_FILE = args.output
    Config.P_VALUE_THRESHOLD = args.pval
    Config.LOGFC_THRESHOLD = args.logfc
    #Config.DAVID_EMAIL = args.email
    #Config.DAVID_API_KEY = args.apikey
    
    # Create and run pipeline
    pipeline = RNASeqAnalysisPipeline(args.disease, args.input)
    pipeline.run_pipeline()

if __name__ == "__main__":
    main()