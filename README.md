# Cancer Subtype Prediction Pipeline (RNA-seq/ Breast cancer)

## Overview

This repository contains an end-to-end **RNA-seq–based cancer subtype prediction pipeline** implemented in **R**. The workflow demonstrates how to process raw count data from GEO, perform normalization and gene annotation, and predict **breast cancer molecular subtypes** using multiple complementary approaches:

* **PAM50 classifier (PAM50 package)**
* **PAM50 via genefu**
* **Signature-based GSVA / ssGSEA scoring**

---

## Key Features

* Automated **GEO data integration** (metadata + raw counts)
* Robust **low-expression gene filtering**
* **edgeR + limma-voom** normalization
* Gene annotation using **EnsDb.Hsapiens.v86**
* Multiple **independent subtype prediction strategies**
* Publication-ready **heatmap visualizations**

---

## Workflow Summary

1. **Data Acquisition**

   * Metadata fetched directly from GEO using `GEOquery`
   * Raw counts loaded from supplementary files

2. **Preprocessing & Normalization**

   * Removal of duplicated genes
   * Filtering genes expressed in ≥80% of samples
   * TMM normalization (`edgeR`)
   * Variance modeling with `voom`

3. **Gene Annotation**

   * Mapping gene symbols to Entrez IDs using `EnsDb.Hsapiens.v86`
   * Harmonization of expression and annotation tables

4. **Subtype Prediction Methods**

   ### Method 1: PAM50 (PAM50 package)

   * Direct application of the PAM50 centroid-based classifier
   * Outputs subtype calls and expression heatmaps

   ### Method 2: PAM50 (genefu)

   * Uses `molecular.subtyping()` with robust PAM50 models
   * Provides subtype probabilities and class assignments

   ### Method 3: GSVA / ssGSEA

   * Signature-based scoring at the single-sample level
   * Flexible framework for custom subtype or pathway gene sets

5. **Visualization & Output**

   * Heatmaps for subtype probabilities and GSVA scores
   * CSV files containing subtype predictions

---

## Dependencies

### R Packages

```
GEOquery
edgeR
limma
data.table
dplyr
caret
genefu
AnnotationDbi
EnsDb.Hsapiens.v86
GSVA
pheatmap
PAM50

