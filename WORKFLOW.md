# postgwas-tools — Workflows

This document provides hands-on examples of how to use the tools included in the `postgwas-tools` package.

---

## Multi-Model Manhattan Plot

Generate a combined Manhattan plot for multiple GWAS result files:

```bash
manhattan_plot -p /data/gwas1.sumstats /data/gwas2.sumstats /data/gwas3.sumstats \
               -t manhattan \
               -o ./results
```

or using a wildcard:
```bash
manhattan_plot -p "/data/*.sumstats" \
               -t manhattan \
               -o ./results
```

**Example output:**

![Manhattan Example](docs/images/manhattan_example.png)

---

## QQ Plot

GWAS inflation
Inspired from (https://github.com/ShujiaHuang/qmplot):

```bash
QQ_plot -p /data/gwas1.sumstats \
        -o ./results
```

![QQ Example](docs/images/qq_example.png)

---

## Identify Lead SNPs (FUMA like, with PLINK1.9)

```bash
leadSNP --sumstats V1_32_SCsylv_l.sumstats \
        --col-a1 A1 --col-a2 A2 \
        --bfile filt_imputed_autosomes_maf-0.01 \
        --plink $PLINK \
        --out tests/V1_32_SCsylv_l
```

This command will extract independent lead SNPs using FUMA-like logic.

---

## Replication Check by rsID

```bash
replication_by_rsid -p /data/discovery.sumstats /data/replication.sumstats -o ./replication_results
```

---

## LDSC Results Extraction

```bash
extract_gencorr -i /data/ldsc/genetic_correlations.log -o ./results
extract_h2 -i /data/ldsc/heritability.log -o ./results
```

---

## Notes

- All plots are saved as `.png` files with **transparent backgrounds** for easy embedding.  
- Works with both `.sumstats` and compressed `.gz` files.  
- Regex paths and file lists are both supported.

---

## Back to Overview

Return to the main [README.md](README.md).