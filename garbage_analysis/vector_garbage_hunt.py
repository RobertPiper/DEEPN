#!/usr/bin/env python3
"""Two complementary, DESeq2-stats-based approaches to finding "vector
alone" garbage genes (non-specific/promiscuous binders), replacing the old
HuORF Garbage.docx list's crude 2.5-fold heuristic:

(A) FREE: pull the Vector-Selected-vs-Non-Selected DESeq2 contrast already
    computed once per Rab as part of every SM6 file's normal collation (36
    independent tests, no new computation) - a gene consistently
    significant there, across many different Rabs' own vector pairs, is
    suspicious regardless of any bait.

(B) LOO ("pseudo-bait"): there are 16 distinct Vector Selected/Non-Selected
    replicate pairs shared across the 36 Rabs (reused 2-9x each). For each
    of the 16, treat it as a held-out "pseudo-bait" and run the existing
    Bait-vs-Vector DESeq2 contrast against the POOLED other 15 pairs as
    the vector group (the R script's infer_bait_name already supports an
    arbitrary number of semicolon-joined replicate columns - no R changes
    needed). A gene significant across most/all 16 held-out iterations,
    independent of which specific replicate is left out, is strong
    evidence of a generic sticky binder rather than any real interaction.

Union of (A) and (B) = the expanded garbage candidate list. MCMC/Bayesian
skipped entirely per instruction - DESeq2 only.
"""
import sys, os, csv, time
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/Users/robertpiper2/LOCAL_DEEPN')
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
os.chdir('/Users/robertpiper2/LOCAL_DEEPN')

import stat_maker_gui_v7 as smg
import sm6_batch as batch
import functions.stat_collation as sc

# the existing 30-gene HuORF Garbage list, for cross-reference only
GARBAGE_GENES = {
    'CYB5R4', 'HSPA8', 'EFNB1', 'SEZ6L', 'ACLY', 'MPDZ', 'STAU2', 'RNF40',
    'CDC42EP1', 'SUGT1', 'KRT222', 'RLIM', 'GOPC', 'CRX', 'PAX7', 'RNF2',
    'POLK', 'ZNF846', 'VCAM1', 'TXNL1', 'TBCCD1', 'BCORL1', 'PCHDGB4',
    'ZNF561', 'ERGIC3', 'RWDD1', 'LRRN1', 'EIF2S2', 'POLR2C', 'C10orf2',
}

OUT_DIR = '/Users/robertpiper2/endosomeLAB Dropbox/Rob Piper/SM6_batch_output'
WORK_DIR = '/private/tmp/claude-503/-Users-robertpiper2-Desktop/d0bd16a2-ab29-49c2-9981-9811e0c1149f/scratchpad/vector_garbage_work'
LOG_PATH = os.path.join(WORK_DIR, 'log.txt')
GC = batch.GC
P_THR = 0.01


def log(msg):
    line = "[%s] %s" % (time.strftime('%H:%M:%S'), msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a') as f:
        f.write(line + '\n')


def get_distinct_vector_pairs():
    pairs = {}
    for fn in batch.CANONICAL:
        cfg = batch.get_config(os.path.join(batch.SM, fn))
        for sk, nk in (('Vector_Selected_1', 'Vector_Non-Selected_1'),
                       ('Vector_Selected_2', 'Vector_Non-Selected_2')):
            pairs[(cfg[sk], cfg[nk])] = True
    return sorted(pairs.keys())


# ---------- Approach A: free, already-computed Vector-vs-Vector per Rab ----------
def approach_a():
    log("=== Approach A: per-Rab Vector-Selected-vs-Non-Selected (already computed) ===")
    sig_count = defaultdict(int)
    tested_count = defaultdict(int)
    n_datasets = 0
    for fn in batch.CANONICAL:
        rid = batch.rab_id(fn)
        path = os.path.join(OUT_DIR, 'SM6_%s_23Aug2026.csv' % rid)
        parsed = smg._parse_sm_csv(path)
        n_datasets += 1
        for g, srow in parsed['stats'].items():
            pv = srow.get('pvalue_vector')
            fc = srow.get('log2FoldChange_vector')
            if pv in (None, '', 'NA') or fc in (None, '', 'NA'):
                continue
            tested_count[g] += 1
            if float(pv) <= P_THR and float(fc) > 0:
                sig_count[g] += 1
    log("  %d datasets tested" % n_datasets)
    return sig_count, tested_count, n_datasets


# ---------- Approach B: leave-one-out pseudo-bait ----------
def build_loo_collated(pseudo_bait_counts, pool_counts_list, output_dir, collated_path):
    """pseudo_bait_counts: {'S':dict, 'N':dict} for the held-out pair.
    pool_counts_list: list of {'S':dict,'N':dict} for the other 15 pairs."""
    genes = sorted(pseudo_bait_counts['S'].keys())
    n_pool = len(pool_counts_list)
    with open(collated_path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['data_column', 'note'])
        w.writerow(['raw_count', 'LOO pseudo-bait vs pooled-other-vectors garbage hunt'])
        w.writerow([])
        w.writerow(['CONFIG_START'])
        w.writerow(['enrich_fold_change', 'enrich_p_val', 'normalized', 'output_directory',
                    'sample_names_selected', 'sample_names_background'])
        w.writerow(['none', 0.05, 'FALSE', output_dir, 'Bait1S', 'Bait1N'])
        vec_s_names = ';'.join('Vector%dS' % (i + 1) for i in range(n_pool))
        vec_n_names = ';'.join('Vector%dN' % (i + 1) for i in range(n_pool))
        w.writerow(['none', 0.05, 'FALSE', output_dir, vec_s_names, vec_n_names])
        w.writerow(['DATA_START'])
        data_cols = ['GENE', 'Bait1S', 'Bait1N']
        for i in range(n_pool):
            data_cols += ['Vector%dS' % (i + 1), 'Vector%dN' % (i + 1)]
        w.writerow(data_cols)
        for g in genes:
            row = [g, pseudo_bait_counts['S'].get(g, 0), pseudo_bait_counts['N'].get(g, 0)]
            for pc in pool_counts_list:
                row += [pc['S'].get(g, 0), pc['N'].get(g, 0)]
            w.writerow(row)
    return collated_path


def approach_b():
    log("=== Approach B: leave-one-out pseudo-bait vs pooled-other-vectors ===")
    pairs = get_distinct_vector_pairs()
    log("  %d distinct vector pairs" % len(pairs))
    rscript_exe = sc.find_rscript()

    log("  loading raw counts for all %d vector pairs..." % len(pairs))
    pair_counts = []
    for sel_fn, non_fn in pairs:
        s = sc.load_raw_counts_summed(GC + sel_fn)
        n = sc.load_raw_counts_summed(GC + non_fn)
        pair_counts.append({'S': s, 'N': n})

    sig_count = defaultdict(int)
    tested_count = defaultdict(int)
    os.makedirs(WORK_DIR, exist_ok=True)
    for held_out in range(len(pairs)):
        t0 = time.time()
        run_dir = os.path.join(WORK_DIR, 'loo_%d' % held_out)
        os.makedirs(run_dir, exist_ok=True)
        pseudo_bait = pair_counts[held_out]
        pool = [pair_counts[i] for i in range(len(pairs)) if i != held_out]
        collated_path = os.path.join(run_dir, 'collated_input.csv')
        build_loo_collated(pseudo_bait, pool, run_dir, collated_path)
        log("  [%d/%d] running DESeq2 (held out: %s)..." % (held_out + 1, len(pairs), pairs[held_out][0]))
        sc.run_r_script(smg.R_SCRIPT_PATH, collated_path, rscript_exe)
        stats = sc.load_stats_output(os.path.join(run_dir, 'everything_combined.csv'))
        n_tested_this_run = 0
        n_sig_this_run = 0
        for g, srow in stats.items():
            pv = srow.get('pvalue_bait1_vs_vector')
            fc = srow.get('log2FoldChange_bait1_vs_vector')
            if pv in (None, '', 'NA') or fc in (None, '', 'NA'):
                continue
            tested_count[g] += 1
            n_tested_this_run += 1
            if float(pv) <= P_THR and float(fc) > 0:
                sig_count[g] += 1
                n_sig_this_run += 1
        log("  [%d/%d] done in %.0fs - %d genes tested, %d significant"
            % (held_out + 1, len(pairs), time.time() - t0, n_tested_this_run, n_sig_this_run))
    return sig_count, tested_count, len(pairs)


def main():
    os.makedirs(WORK_DIR, exist_ok=True)
    open(LOG_PATH, 'w').close()

    sig_a, tested_a, n_a = approach_a()
    sig_b, tested_b, n_b = approach_b()

    KRT_check_a = sorted((g, sig_a[g], tested_a.get(g, 0)) for g in sig_a if g.startswith('KRT'))
    KRT_check_b = sorted((g, sig_b[g], tested_b.get(g, 0)) for g in sig_b if g.startswith('KRT'))
    log("\nKRT** genes, Approach A (of %d Rab datasets): %s" % (n_a, KRT_check_a))
    log("KRT** genes, Approach B (of %d LOO iterations): %s" % (n_b, KRT_check_b))

    # Three tracked sets, exactly as specified: the union from Approach A
    # alone, the union from Approach B alone, and their combined union -
    # kept as separate, explicit columns rather than merged into one
    # opaque list.
    set_a = set(sig_a.keys())
    set_b = set(sig_b.keys())
    set_union = set_a | set_b
    log("\nApproach A union: %d genes | Approach B union: %d genes | Combined union: %d genes"
        % (len(set_a), len(set_b), len(set_union)))

    out_csv = os.path.join(OUT_DIR, 'Rab_vector_garbage_candidates.csv')
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['gene', 'in_approachA_union', 'approachA_sig_of_%d_Rabs' % n_a, 'approachA_frac',
                    'in_approachB_union', 'approachB_sig_of_%d_LOO' % n_b, 'approachB_frac',
                    'in_combined_union', 'on_old_garbage_list'])
        for g in sorted(set_union):
            a_n = sig_a.get(g, 0); a_t = tested_a.get(g, 0)
            b_n = sig_b.get(g, 0); b_t = tested_b.get(g, 0)
            w.writerow([g, g in set_a, a_n, round(a_n / a_t, 3) if a_t else '',
                        g in set_b, b_n, round(b_n / b_t, 3) if b_t else '',
                        g in set_union, g in GARBAGE_GENES])
    log("Wrote %s" % out_csv)
    log("=== DONE ===")


if __name__ == '__main__':
    main()
