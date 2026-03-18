# RNA-seq Analysis Pipeline (Automated)

## Overview
This project was developed as part of a class assignment where we were asked to perform an RNA-seq analysis pipeline manually — starting from identifying differentially expressed genes to exploring pathways and drug targets. Our teacher also suggested that this workflow could be automated using a script.
Instead of following the traditional manual approach, my friend and I decided to implement the entire pipeline programmatically. As bioinformatics students, we aimed to build a solution that is faster, reproducible, and scalable.

## What This Pipeline Does
This Python-based pipeline automates the complete RNA-seq downstream analysis workflow:

* Takes a **TSV file (e.g., GEO2R output)** as input
* Filters **Differentially Expressed Genes (DEGs)** using statistical thresholds
* Performs **pathway enrichment analysis** using Enrichr
* Identifies **disease-associated genes**
* Maps **drug–gene interactions**
* Generates a structured **Excel report** with all results

## Workflow

```
RNA-seq input  
   ↓  
DEG filtering  
   ↓  
Pathway enrichment (Enrichr)  
   ↓  
Disease genes  
   ↓  
Drug targets  
   ↓  
Excel report  
```

## Why This Approach
Manual analysis of RNA-seq data involves multiple tools, repeated steps, and a high chance of human error. This pipeline addresses those issues by:
* Automating the entire workflow
* Reducing analysis time from hours to minutes
* Ensuring reproducibility
* Making it easy to apply the same analysis to any disease dataset

## Key Features
* Dynamic detection of input columns (gene, p-value, logFC)
* Configurable thresholds for DEG filtering
* Integration with biological databases (Enrichr, KEGG)
* Modular functions for easy extension and customization
* Clean Excel output with multiple result sheets

## Use Case
Although demonstrated using **breast cancer data**, this pipeline is **not disease-specific**. It can be used for any condition by simply providing a relevant RNA-seq dataset and disease name.

**How to Run**
python venom.py "breast cancer" --input input.tsv
The code file and the input .tsv file should be within the same folder the run the above mentioned command in command Prompt(cmd)

**Future Perspective**
The pipeline can be improved by integrating more comprehensive and updated databases such as Reactome, WikiPathways, or DrugBank for better pathway and drug analysis. Additionally, using adaptive thresholds and combining multiple data sources could increase accuracy and biological relevance.

---
