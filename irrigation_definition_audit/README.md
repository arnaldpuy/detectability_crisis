# Irrigation Definition Audit Pipeline

This repository provides a reproducible pipeline to extract and audit how global irrigation datasets define irrigation across scientific studies.

## Overview

Global irrigation maps differ not only in magnitude, but in how irrigation is defined. This repository implements a structured audit of these definitions using a large language model (DeepSeek-R1, 32B parameters).

Each document is classified across six dimensions:

* **Irrigation target** (e.g., area equipped vs actually irrigated)
* **Representation unit** (binary presence vs area magnitude)
* **Temporal concept** (annual, seasonal, snapshot)
* **Treatment of paddy rice**
* **Treatment of permanent crops**
* **Irrigation basis** (infrastructure-based vs water-use based)

The pipeline enforces strict, evidence-based classification and supports reproducible analysis of definitional heterogeneity across datasets.

---

## Methodological features

* **Constrained classification**: labels must be selected from predefined categories
* **Evidence-based extraction**: all classifications require explicit textual support
* **Conservative inference**: absence of evidence results in `"Unclear"`
* **Adversarial verification**: a second model pass validates each classification
* **Retrieval auditing**: relevance of extracted pages is quantified (coverage and recall proxy)

---

## Repository structure

```
.
├── src/
│   └── extract_irrigation_definitions.py
├── prompts/
│   └── prompt.txt
├── outputs/
│   └── irrigation_definition_table.csv
├── requirements.txt
├── README.md
```

---

## Installation

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Usage

Run the extraction pipeline:

```bash
python src/extract_irrigation_definitions.py \
  --pdf_dir ./pdfs \
  --cache_dir ./cache \
  --out_dir ./outputs
```

To prevent system sleep during long runs (macOS):

```bash
caffeinate -dimsu python src/extract_irrigation_definitions.py ...
```

---

## Outputs

The pipeline produces:

* `irrigation_definition_table.csv` → main analysis dataset
* `irrigation_definition_table.xlsx` → spreadsheet version
* `irrigation_definition_details.jsonl` → full structured outputs (including citations)

Each record includes:

* assigned label per category
* confidence score
* supporting quote and page reference
* retrieval quality metrics

---

## Prompt transparency

The full prompt used for classification is provided in:

```
prompts/prompt.txt
```

The prompt enforces strict label selection, explicit evidence requirements, and conservative classification rules.

---

## Reproducibility

* Deterministic outputs (temperature = 0)
* Fully specified prompt and schema
* Cached intermediate results
* Retrieval and classification steps are auditable

---

## Data availability

The extracted irrigation definition dataset is provided in this repository as a CSV file.

Input PDFs are not distributed due to licensing constraints. Users should provide their own corpus of documents.

---

## Citation

If you use this pipeline, please cite:

Puy, A. et al. (2026). Code of Where irrigation exists is globally contested. Zenodo. doi: 10.5281/zenodo.19001232.

---

## License

MIT License

Copyright (c) 2026 Arnald Puy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
