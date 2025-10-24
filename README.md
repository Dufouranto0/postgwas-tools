#  postgwas-tools

**postgwas-tools** is a Python package that provides a collection of tools for **post-GWAS (Genome-Wide Association Study) analysis**, annotation, and visualization.  
It simplifies steps such as plotting, SNP annotation, and result parsing for large-scale GWAS projects.

---

## Installation

### Option 1 — From source
```bash
git clone https://github.com/Dufouranto0/postgwas-tools.git
cd postgwas-tools
python3 -m venv venv
. venv/bin/activate
pip install -e .
```

### Option 2 — From PyPI (planned)
```bash
pip install postgwas-tools
```

---

## Package Structure

```
src/postgwas_tools/
├── annot/
│   ├── fuma/                 # FUMA-related post-GWAS visualization
│   ├── ldsc/                 # LD Score Regression utilities
│   ├── magma/                # MAGMA annotation helpers
│   ├── replication/          # Replication and cross-cohort checks
│   └── utils.py              # Common helper functions
└── mostest/                  # MOSTest-specific result processors
```

---

## Quick Start Example

You can generate a multi-model Manhattan plot directly from the command line:

```bash
manhattan_plot -p /path/to/gwas1.sumstats /path/to/gwas2.sumstats \
               -t manhattan \
               -o ./results
```

Or using a file pattern:
```bash
manhattan_plot -p "/path/to/*.sumstats"
```

This will create a high-resolution `manhattan_plot.png` in the output directory.

If you just want a light manhattan plot, you can do:
```bash
light_manhattan_plot -p /path/to/gwas_input.txt
```

---

## Documentation & Workflows

For detailed usage examples, including **QQ plots**, **locus-based analysis**, and **replication workflows**, see:
 [WORKFLOW.md](WORKFLOW.md)

---

## Dependencies

- Python ≥ 3.9  
- pandas, numpy, matplotlib, seaborn, scipy, qmplot, h5py

All dependencies are automatically installed when you run `pip install -e .`.

---

## Citation

If you use `postgwas-tools` in your work, please cite:

> Dufournet A. *postgwas-tools: a Python suite for streamlined post-GWAS analysis*. CEA NeuroSpin, 2025.

---

## Links

- [LICENSE](LICENSE)  
- [Homepage](https://github.com/Dufouranto0/postgwas-tools)  
- [Issue Tracker](https://github.com/Dufouranto0/postgwas-tools/issues)