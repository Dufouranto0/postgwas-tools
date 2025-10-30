#!/usr/bin/env python3
"""
summarize_replication.py

Aggregate region-wise replication results.

Usage:
    python3 summarize_replication.py \
        --discovery /path/to/UKB/Champollion_V1_32 \
        --replication /path/to/ABCD/Champollion_V1_32 \
        --out /path/to/final_replication_summary.tsv \
        --log /path/to/replication_summary.log

It will:
- Loop over all regions that contain a replication_matches.txt file
- Merge discovery leadSNPs.txt and replication results by rsID
- Compute replication success rate (Bonferroni-corrected 0.05/N)
- Output a region-level summary TSV and a global summary log
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def summarize_region(lead_file, repl_file):
    """Return merged DataFrame and replication statistics."""
    lead = pd.read_csv(lead_file, sep="\t", dtype=str)
    repl = pd.read_csv(repl_file, sep="\t", dtype=str)

    rs_col_lead = next((c for c in lead.columns if c.lower() in ("rsid","rs_id","snp","variant","markername")), None)
    rs_col_repl = next((c for c in repl.columns if c.lower() in ("rsid","rs_id","snp","variant","markername")), None)
    p_col_repl = next((c for c in repl.columns if c.lower() in ("pval","p_value","p-value","p")), None)

    if not rs_col_lead or not rs_col_repl or not p_col_repl:
        raise ValueError(f"Missing rsID or PVAL columns in {lead_file} or {repl_file}")

    cols_lead_keep = [c for c in lead.columns if c.lower() in ("uniqid","rsid","chr","pos","p","a1","a2")]
    lead = lead[cols_lead_keep].rename(columns={"p": "p_dis"})

    cols_repl_keep = [rs_col_repl, p_col_repl, "A1", "A2"]
    cols_repl_keep = [c for c in cols_repl_keep if c in repl.columns]
    repl = repl[cols_repl_keep].rename(columns={p_col_repl: "p_rep"})

    merged = pd.merge(lead, repl, left_on=rs_col_lead, right_on=rs_col_repl, how="inner")
    n_tested = len(merged)
    if n_tested == 0:
        return merged, 0, 0, 0, 0

    merged["p_rep"] = pd.to_numeric(merged["p_rep"], errors="coerce")
    threshold = 0.05 / n_tested
    merged["replicated"] = merged["p_rep"] < threshold

    n_rep = merged["replicated"].sum()
    perc = (n_rep / n_tested) * 100
    return merged, n_tested, n_rep, threshold, perc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--discovery", required=True)
    parser.add_argument("--replication", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    results = {"White": [], "All": []}

    with open(args.log, "w") as log:
        log.write("Population\tRegion\tN_tested\tN_replicated\tThreshold\tReplication_%\n")

        for root, _, files in os.walk(args.replication):
            if "replication_matches.txt" not in files:
                continue

            repl_file = os.path.join(root, "replication_matches.txt")
            region_path = "/".join(os.path.relpath(repl_file, args.replication).split("/")[:2])

            # Try detecting population
            if "White" in root or "white" in root:
                pop = "White"
            else:
                pop = "All"

            # Derive matching discovery file (may differ by structure)
            lead_file = os.path.join(args.discovery, region_path, "NOPCA/white.British.ancestry/FUMA/leadSNPs.txt")
            if not os.path.exists(lead_file):
                print(f"Missing discovery file for {region_path}")
                continue

            try:
                merged, n_tested, n_rep, thr, perc = summarize_region(lead_file, repl_file)
                merged["region"] = region_path
                merged["population"] = pop
                results[pop].append({
                    "region": region_path,
                    "N_tested": n_tested,
                    "N_replicated": n_rep,
                    "Threshold": thr,
                    "Replication_%": perc
                })
                log.write(f"{pop}\t{region_path}\t{n_tested}\t{n_rep}\t{thr:.2e}\t{perc:.2f}\n")
                print(f"{pop} | {region_path}: {n_rep}/{n_tested} ({perc:.1f}%)")

            except Exception as e:
                print(f"{region_path} failed: {e}")

    # Save per-population results
    for pop, data in results.items():
        if not data:
            continue
        df = pd.DataFrame(data)
        pop_out = args.out.replace(".tsv", f"_{pop.lower()}.tsv")
        df.to_csv(pop_out, sep="\t", index=False)
        print(f"Saved {pop} summary: {pop_out}")

    # Create comparison plot if both populations available
    if all(len(results[p]) > 0 for p in ["White", "All"]):
        df_white = pd.DataFrame(results["White"])
        df_all = pd.DataFrame(results["All"])
        merged_plot = pd.merge(df_white, df_all, on="region", suffixes=("_White", "_All"))
        merged_plot["region"] = merged_plot["region"].apply(lambda x : x.split('/name')[0])
        merged_plot = merged_plot.sort_values(by="region")

        plt.figure(figsize=(10, 6))
        x = range(len(merged_plot))
        plt.bar([i - 0.2 for i in x], merged_plot["Replication_%_White"], width=0.4, label="White")
        plt.bar([i + 0.2 for i in x], merged_plot["Replication_%_All"], width=0.4, label="All")

        plt.xticks(x, merged_plot["region"], rotation=80, ha="right", fontsize=8)
        plt.ylabel("Replication %")
        plt.title("Replication Success by Region and Population")
        plt.legend()
        plt.tight_layout()

        plot_path = args.out.replace(".tsv", "_comparison_barplot.png")
        plt.savefig(plot_path, dpi=300)
        print(f"Saved bar plot: {plot_path}")


if __name__ == "__main__":
    main()
