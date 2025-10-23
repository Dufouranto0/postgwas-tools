# coding: utf-8
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# Antoine Dufournet
##########################################################################
"""
leadSNP.py

Usage example:
python3 leadSNP.py \
  --sumstats most_orig.sumstats \
  --col-snp SNP --col-chr CHR --col-bp BP --col-a1 A1 --col-a2 A2 --col-p PVAL \
  --bfile /ccc/scratch/cont003/n4h00001/dufourna/25irene_AD_UKB_TIV/Champollion_V1_32/White/gen/filt_imputed_autosomes_maf-0.01 \
  --plink "pcocc-rs run n4h00001rs:plink1.9 plink1.9 --" \
  --out results

Defaults:
 - clump r2 thresholds: 0.1 (loci), 0.6 (independent)
 - clump-kb: 1000 (kb)
 - merge-distance: 250 (kb)
"""

import argparse
import os
import subprocess
import sys
import pandas as pd
import gc
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description="Run PLINK clumping and produce a FUMA-like leadSNPs table.")
    p.add_argument("--sumstats", required=True, help="GWAS summary stats file (tab/space-delimited, may have header).")
    p.add_argument("--col-snp", default="SNP", help="Column name for SNP id in sumstats (default: SNP)")
    p.add_argument("--col-chr", default="CHR", help="Column name for chromosome (default: CHR)")
    p.add_argument("--col-bp", default="BP", help="Column name for base pair position (default: BP)")
    p.add_argument("--col-a1", default="A1", help="Column name for effect allele / A1 (default: A1)")
    p.add_argument("--col-a2", default="A2", help="Column name for other allele / A2 (default: A2)")
    p.add_argument("--col-p", default="P", help="Column name for p-value (default: P)")
    p.add_argument("--bfile", required=True, help="PLINK --bfile prefix (path to .bed/.bim/.fam without extension)")
    p.add_argument("--plink", default="plink", help="PLINK binary to call (default 'plink')")
    p.add_argument("--r2-loci", type=float, default=0.1, help="r2 for loci clumping (default: 0.1)")
    p.add_argument("--r2-ind", type=float, default=0.6, help="r2 for independent SNPs clumping (default: 0.6)")
    p.add_argument("--clump-p1", type=float, default=5e-8, help="PLINK --clump-p1 (default 5e-8)")
    p.add_argument("--clump-p2", type=float, default=5e-6, help="PLINK --clump-p2 (default 5e-6)")
    p.add_argument("--clump-kb", type=int, default=1000, help="PLINK --clump-kb (default 1000)")
    p.add_argument("--merge-distance", type=int, default=250, help="distance in kb to merge loci (default 250)")
    p.add_argument("--out", default="results", help="Output directory (default 'results')")
    p.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate reformatted sumstats file")
    return p.parse_args()

def read_sumstats(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        df = pd.read_csv(path, sep="\s+")
    return df

def write_good_format(df, mapping, outpath, p2_threshold):
    # ensure columns exist
    for c in mapping.values():
        if c not in df.columns:
            raise KeyError(f"Column '{c}' not found in sumstats input. Available columns: {list(df.columns)}")
    # Build minimal table with needed columns
    sumstats = pd.DataFrame({
        "SNP": df[mapping["snp"]],
        "CHR": df[mapping["chr"]],
        "BP": df[mapping["bp"]],
        "A1": df[mapping["a1"]],
        "A2": df[mapping["a2"]],
        "P": df[mapping["p"]]
    })

    # drop NA in essential fields
    sumstats = sumstats.dropna(subset=["SNP", "CHR", "BP", "P"])

    # Apply p2 threshold prefilter
    if p2_threshold is not None:
        sumstats = sumstats[sumstats["P"] <= float(p2_threshold)].copy()

    # Create uniqID here with alphabetically sorted alleles
    def make_uid(row):
        a1, a2 = sorted([str(row["A1"]), str(row["A2"])])
        return f"{int(row['CHR'])}:{int(row['BP'])}:{a1}:{a2}"
    sumstats["uniqID"] = sumstats.apply(make_uid, axis=1)

    # Save and return
    sumstats.to_csv(outpath, sep="\t", index=False)
    return sumstats


def run_plink_clump(plink_bin, bfile_prefix, sumstats_file, out_prefix, p1, p2, r2, kb):
    cmd = [
        plink_bin,
        "--bfile", bfile_prefix,
        "--clump", sumstats_file,
        "--clump-p1", str(p1),
        "--clump-p2", str(p2),
        "--clump-r2", str(r2),
        "--clump-kb", str(kb),
        "--out", out_prefix
    ]
    print("Running PLINK:", " ".join(cmd))
    res = subprocess.run(" ".join(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    if res.returncode != 0:
        print("PLINK failed. stdout:\n", res.stdout)
        print("PLINK stderr:\n", res.stderr, file=sys.stderr)
        raise RuntimeError("PLINK clumping failed.")
    else:
        # optionally print trimmed stdout
        print(res.stdout.splitlines()[-5:] if res.stdout else "PLINK completed.")
    return out_prefix + ".clumped"

def read_clumped(path):
    # PLINK clumped files typically have whitespace-separated columns and a header row.
    if not os.path.exists(path):
        print(f"Warning: clumped file {path} not found; returning empty DataFrame")
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\s+", comment="#")
    except Exception:
        print("Issue when reading plink output")
    return df

def merge_loci_by_distance(loci_df, merge_dist):
    """Merge PLINK loci by distance threshold. If merge_dist == 0, return loci as-is."""
    if loci_df.empty:
        return [], loci_df

    loci_df = loci_df.copy()
    loci_df['CHR'] = loci_df['CHR'].astype(str)
    loci_df['BP'] = loci_df['BP'].astype(int)
    loci_df['P'] = loci_df['P'].astype(float)
    loci_df = loci_df.sort_values(['CHR', 'BP']).reset_index(drop=True)

    if merge_dist <= 0:
        # Keep each PLINK clump locus as-is — no merging
        merged = []
        for i, r in loci_df.iterrows():
            merged.append({
                'chr': r['CHR'],
                'start': r['BP'],
                'end': r['BP'],
                'indices': [i],
                'lead_row': r,
                'lead_local_idx': i,
                'subset_df': pd.DataFrame([r])
            })
        return merged, loci_df

    # --- Otherwise merge loci by distance ---
    groups = []
    current = {'chr': loci_df.loc[0,'CHR'], 'start': loci_df.loc[0,'BP'],
               'end': loci_df.loc[0,'BP'], 'indices': [0]}
    for i in range(1, len(loci_df)):
        r = loci_df.loc[i]
        if r['CHR'] == current['chr'] and r['BP'] - current['end'] <= merge_dist:
            current['indices'].append(i)
            current['end'] = max(current['end'], r['BP'])
        else:
            groups.append(current)
            current = {'chr': r['CHR'], 'start': r['BP'], 'end': r['BP'], 'indices':[i]}
    groups.append(current)

    merged = []
    for g in groups:
        subset = loci_df.loc[g['indices']]
        lead_local_idx = subset['P'].idxmin()
        lead_row = subset.loc[lead_local_idx]
        merged.append({
            'chr': g['chr'],
            'start': g['start'],
            'end': g['end'],
            'indices': g['indices'],
            'lead_row': lead_row,
            'lead_local_idx': lead_local_idx,
            'subset_df': subset
        })
    return merged, loci_df


def sort_by_chr_pos(df):
    """
    Sort numerically by chromosome and position extracted from uniqID
    """
    df[['chr_sort', 'pos_sort']] = df['uniqID'].str.extract(r'^(\d+):(\d+):')
    df['chr_sort'] = df['chr_sort'].astype(int)
    df['pos_sort'] = df['pos_sort'].astype(int)
    df = df.sort_values(['chr_sort', 'pos_sort']).reset_index(drop=True)
    df = df.drop(columns=['chr_sort', 'pos_sort'])
    return df


def parse_SP2_field(sp2_field):
    """Parse PLINK SP2 field like 'rs1(1),rs2(1),...' -> ['rs1','rs2',...]"""
    if pd.isna(sp2_field):
        return []
    # split on comma, strip whitespace, extract part before '('
    items = [x.strip() for x in str(sp2_field).split(',') if x.strip()]
    snps = []
    for it in items:
        # Some items might be like 'rs123(1)' or 'NONE'
        if it.upper() == "NONE":
            continue
        if '(' in it:
            snp = it.split('(')[0].strip()
        else:
            snp = it
        if snp:
            snps.append(snp)
    return snps

def build_fuma_tables(merged_loci, loci_df, ind_df, good_sumstats_df, outdir):
    """
    Build both leadSNPs and IndSigSNPs tables.
    Optimized: use SP2 from loci_df and filter SP2 entries to those present in ind_df (r2=0.6 SNPs).
    """
    merged, loci_df_sorted = merged_loci

    # allele_map for uniqID creation and p-values (from good_sumstats_df)
    allele_map = {}
    for _, r in good_sumstats_df.iterrows():
        allele_map[str(r['SNP'])] = (
            str(r['A1']), str(r['A2']), str(r['CHR']),
            int(r['BP']), float(r['P']), str(r['uniqID'])
        )

    # Prepare a fast lookup set of independent SNPs (from PLINK r2=0.6 clump)
    indep_snp_set = set()
    ind_by_snp = {}
    if not ind_df.empty:
        for _, r in ind_df.iterrows():
            s = str(r['SNP'])
            indep_snp_set.add(s)
            ind_by_snp[s] = {'CHR': str(r['CHR']), 'BP': int(r['BP']), 'P': float(r.get('P', r.get('p', float('nan'))))}

    lead_rows = []
    # Build lead SNP (locus) table using SP2 from PLINK clump_loci
    for m in merged:
        subset = m['subset_df']
        # PLINK clump_loci rows include SP2 field; subset likely has one row per original loci (if not merged)
        # We'll build the set of significant SNPs as: lead + those listed in SP2 that also appear in indep_snp_set
        lead_snp = str(m['lead_row']['SNP'])
        lead_p = float(m['lead_row']['P'])
        lead_chr = str(m['lead_row']['CHR'])
        lead_bp = int(m['lead_row']['BP'])

        inds = [lead_snp]

        # For each row in the subset, parse its SP2 field and collect candidate SNPs
        # (This handles the case where a merged locus contains multiple original loci each with SP2)
        candidate_sp2 = []
        for _, row in subset.iterrows():
            sp2_field = row.get('SP2') if 'SP2' in row.index else None
            candidate_sp2.extend(parse_SP2_field(sp2_field))

        # Keep only those SP2 SNPs that are present in the independent r2=0.6 clump result
        for s in candidate_sp2:
            if s == lead_snp:
                continue
            if s in indep_snp_set:
                inds.append(s)

        # uniqID: use allele_map if available, with alphabetical allele order
        a1, a2 = ("NA", "NA")
        if lead_snp in allele_map:
            a1_raw, a2_raw, chr_, bp_, _, _ = allele_map[lead_snp]
            a1, a2 = sorted([a1_raw, a2_raw])
        uniqid = allele_map[lead_snp][5] if lead_snp in allele_map else f"{lead_chr}:{lead_bp}:{a1}:{a2}"

        # Count total SP2 SNPs (including the lead) from the PLINK locus table
        sp2_all = []
        for _, row in subset.iterrows():
            sp2_all.extend(parse_SP2_field(row.get("SP2")))
        sp2_count = len(set([lead_snp] + sp2_all))

        lead_rows.append({
            "chr": lead_chr,
            "pos": lead_bp,
            "rsID": lead_snp,
            "p": lead_p,
            "uniqID": uniqid,
            "IndSigSNPs": ";".join(inds),
            "nIndSigSNPs": len(inds),
            "nSP2SNPs": sp2_count  # <-- new field tracking number of SP2 SNPs
        })

    # free the big good_sumstats_df if you won't need it anymore for memory
    # (we used allele_map to capture what we needed)
    try:
        del good_sumstats_df
        gc.collect()
    except Exception:
        pass

    # Sort loci by chr:pos and number them
    lead_df = pd.DataFrame(lead_rows)
    lead_df = sort_by_chr_pos(lead_df)
    lead_df.insert(0, "No", range(1, len(lead_df) + 1))
    lead_df.insert(1, "GenomicLocus", range(1, len(lead_df) + 1))
    lead_df = lead_df[["No","GenomicLocus","uniqID","rsID","chr","pos","p","nIndSigSNPs","IndSigSNPs"]]

    lead_path = os.path.join(outdir, "leadSNPs.txt")
    lead_df.to_csv(lead_path, sep="\t", index=False)
    print(f"Lead SNP table saved to: {lead_path}")

    # ----------------------------
    # Build independent SNP table from lead_df IndSigSNPs lists
    # ----------------------------
    ind_rows = []
    # For quick retrieval of allele/p/pos from the allele_map we built earlier, rebuild local small mapping
    # Note: allele_map was built above; if it's big you can instead read per-snp from disk or other source
    for _, lead in lead_df.iterrows():
        locus_id = int(lead['GenomicLocus'])
        inds = [s for s in str(lead['IndSigSNPs']).split(";") if s]
        sp2_count = int(lead.get("nSP2SNPs", len(inds)))
        for s in inds:
            # prefer allele info from the allele_map; if missing, fallback to ind_by_snp
            if s in allele_map:
                a1_raw, a2_raw, chr_, bp_, pval, uniqid = allele_map[s]
                a1, a2 = sorted([a1_raw, a2_raw])
            elif s in ind_by_snp:
                info = ind_by_snp[s]
                chr_ = info['CHR']
                bp_ = info['BP']
                pval = info['P']
                uniqid = f"{chr_}:{bp_}:NA:NA"
            else:
                chr_, bp_, pval = ("NA", "NA", float('nan'))
                uniqid = f"{chr_}:{bp_}:NA:NA"
            ind_rows.append({
                "No": len(ind_rows) + 1,
                "GenomicLocus": locus_id,
                "uniqID": uniqid,
                "rsID": s,
                "chr": chr_,
                "pos": bp_,
                "p": pval,
                "nSNPs": sp2_count,
                "nGWASSNPs": sp2_count
            })

    ind_df_out = pd.DataFrame(ind_rows)
    ind_df_out = sort_by_chr_pos(ind_df_out)
    ind_path = os.path.join(outdir, "IndSigSNPs.txt")
    ind_df_out.to_csv(ind_path, sep="\t", index=False)
    print(f"Independent SNP table saved to: {ind_path}")

    return lead_df, ind_df_out



def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    df = read_sumstats(args.sumstats)
    mapping = {
        "snp": args.col_snp,
        "chr": args.col_chr,
        "bp": args.col_bp,
        "a1": args.col_a1,
        "a2": args.col_a2,
        "p": args.col_p
    }
    okstats_path = os.path.join(args.out, "good_format_gwas.txt")
    # prefilter by clump_p2 to reduce file size/time/memory
    good_df = write_good_format(df, mapping, okstats_path, p2_threshold=args.clump_p2)
    print(f"Reformatted (prefiltered p<={args.clump_p2}) sumstats written to {okstats_path}; {len(good_df)} rows")
    # free the original large dataframe
    try:
        del df
        gc.collect()
    except Exception:
        pass

    clump_loci_prefix = os.path.join(args.out, "clump_loci")
    clump_ind_prefix = os.path.join(args.out, "clump_ind")

    clumped_loci_path = run_plink_clump(args.plink, 
                                        args.bfile, 
                                        okstats_path, 
                                        clump_loci_prefix, 
                                        args.clump_p1, 
                                        args.clump_p2, 
                                        args.r2_loci, 
                                        args.clump_kb)
    clumped_ind_path = run_plink_clump(args.plink, 
                                        args.bfile, 
                                        okstats_path, 
                                        clump_ind_prefix, 
                                        args.clump_p1, 
                                        args.clump_p2, 
                                        args.r2_ind, 
                                        args.clump_kb)

    loci_df = read_clumped(clumped_loci_path)
    ind_df = read_clumped(clumped_ind_path) 

    for d in (loci_df, ind_df): 
        if not d.empty and 'P' not in d.columns and 'p' in d.columns: 
            d.rename(columns={'p': 'P'}, inplace=True)

    # --- Step 1: generate "no-merge" loci outputs (r2=0.1 lead SNPs, r2=0.6 independent SNPs) ---
    merged_loci_nom = merge_loci_by_distance(loci_df, 0)  # no merging
    lead_df_nom, ind_df_nom = build_fuma_tables(merged_loci_nom, merged_loci_nom[1], ind_df, good_df, args.out) 

    # --- Step 2: create merged GenomicRiskLoci.txt (using merge_distance) ---
    #merged_loci_merged = merge_loci_by_distance(loci_df, args.merge_distance * 1000)
    #merged_groups, merged_df = merged_loci_merged
    #loci_df_out = pd.DataFrame(loci_rows)
    #loci_path = os.path.join(args.out, "GenomicRiskLoci.txt")
    #loci_df_out.to_csv(loci_path, sep="\t", index=False)
    #print(f"GenomicRiskLoci file saved to: {loci_path}")

    # --- Step 3: save parameter file ---
    config_path = os.path.join(args.out, "params.config")
    with open(config_path, "w") as f:
        for k, v in vars(args).items():
            f.write(f"{k}={v}\n")
    print(f"Parameters saved to: {config_path}")

    # --- Step 4: cleanup ---
    if not args.keep_intermediate:
        try:
            os.remove(okstats_path)
        except Exception:
            pass

if __name__ == "__main__":
    main()
