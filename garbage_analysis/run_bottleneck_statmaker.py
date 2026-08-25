#!/usr/bin/env python3
"""Runs the full Stat Maker v7 pipeline (raw+norm DESeq2, Bayesian/MCMC,
junction collation) on the HuORFeome bottleneck pooled sample as Bait1,
paired with its designated Non-Selected partner (pTEFGBD_Vector5_S47_NON -
the bottleneck pool has no true non-selected partner of its own, so this
pairing follows the precedent set in the original 2018
stat_maker_hORF_TEF_control_3way.csv run), against a real Vector1/Vector2
reference pair (Vector11 + Vector2, same as that original run used), no
Bait2. Output: SM6-format rich CSV, same as every other Rab dataset, so
the same _gene_passes_criteria / GUI-default criteria can mine it directly
for bottleneck hits.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/Users/robertpiper2/LOCAL_DEEPN')
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
os.chdir('/Users/robertpiper2/LOCAL_DEEPN')

import stat_maker_gui_v7 as smg
import sm6_batch as batch
import functions.stat_collation as sc

GC = "/Volumes/PiperLabDataDisk/DEEPN_2018_RabGTPase/gene_count_summary/"
BQP_ROOT = "/Volumes/PiperLabDataDisk/DEEPN_2018_RabGTPase"
OUT_DIR = "/Users/robertpiper2/Dropbox (endosomeLAB)/SM6_batch_output"
WORK_DIR = '/private/tmp/claude-503/-Users-robertpiper2-Desktop/d0bd16a2-ab29-49c2-9981-9811e0c1149f/scratchpad/bottleneck_run'
LOG_PATH = os.path.join(WORK_DIR, 'log.txt')

csv_paths = {
    'Bait1S': GC + 'HuORF_pTEF_Control_S202_SEL_summary.csv',
    'Bait1N': GC + 'pTEFGBD_Vector5_S47_NON_summary.csv',
    'Vector1S': GC + 'pTEFGBD_Vector11_S155_SEL_summary.csv',
    'Vector1N': GC + 'pTEFGBD_Vector11_S154_NON_summary.csv',
    'Vector2S': GC + 'pTEFGBD_Vector2_S18_SEL_summary.csv',
    'Vector2N': GC + 'pTEFGBD_Vector2_S17_NON_summary.csv',
}
has_bait2 = False
DATASET_LABELS = {
    'Bait1N': 'Bait1_Non-Selected (Vector5)', 'Bait1S': 'Bait1_Selected (Bottleneck pool)',
    'Vector1N': 'Vector_Non-Selected_1', 'Vector2N': 'Vector_Non-Selected_2',
    'Vector1S': 'Vector_Selected_1', 'Vector2S': 'Vector_Selected_2',
}


def log(msg):
    line = "[%s] %s" % (time.strftime('%H:%M:%S'), msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a') as f:
        f.write(line + '\n')


def main():
    os.makedirs(WORK_DIR, exist_ok=True)
    open(LOG_PATH, 'w').close()

    for k, p in csv_paths.items():
        if not os.path.exists(p):
            raise FileNotFoundError("%s: %s" % (k, p))

    rscript_exe = sc.find_rscript()
    log("=== Bottleneck-vs-Vector Stat Maker v7 run starting ===")

    log("loading raw counts...")
    raw_count_data = {k: sc.load_raw_counts_summed(p) for k, p in csv_paths.items()}
    genes = sorted(raw_count_data['Bait1S'].keys())
    ppm_data = {k: sc.load_ppm_summed(p) for k, p in csv_paths.items()}

    collated_path = os.path.join(WORK_DIR, 'collated_input.csv')
    sc.build_collated_input_v5(raw_count_data, WORK_DIR, collated_path, has_bait2)
    log("running raw-counts + poscounts DESeq2 (p-value_raw)...")
    sc.run_r_script(smg.R_SCRIPT_PATH, collated_path, rscript_exe)
    stats = sc.load_stats_output(os.path.join(WORK_DIR, 'everything_combined.csv'))
    stats = smg._suffix_bait_vs_vector(stats, has_bait2, '_raw')

    log("running normalized-PPM DESeq2 (p-value_norm)...")
    norm_dir = os.path.join(WORK_DIR, 'norm')
    os.makedirs(norm_dir, exist_ok=True)
    collated_norm = os.path.join(norm_dir, 'collated_input.csv')
    sc.build_collated_input_v3(ppm_data, norm_dir, collated_norm, has_bait2)
    sc.run_r_script(smg.R_SCRIPT_PATH, collated_norm, rscript_exe)
    stats_norm = sc.load_stats_output(os.path.join(norm_dir, 'everything_combined.csv'))
    smg._merge_bait_vs_vector_norm(stats, stats_norm, has_bait2)

    log("running Bayesian/MCMC (slow step)...")
    _, bayes_csv, bayes_overdispersion = sc.run_bayesian_r_script(
        smg.BAYESIAN_R_SCRIPT_PATH, csv_paths, smg.BAYESIAN_THRESHOLD_PPM, WORK_DIR, rscript_exe)
    bayesian_stats = sc.load_bayesian_stats_output(bayes_csv)

    log("loading junction (.bqp) data...")
    bqp_data = {}
    for k, p in csv_paths.items():
        basename = os.path.basename(p).replace('_summary.csv', '')
        bqp_path = sc.find_bqp_path(BQP_ROOT, basename)
        bqp_data[k] = sc.load_bqp(bqp_path) if bqp_path else {}
        log("  %s -> %s" % (k, bqp_path))

    out_name = 'SM6_BottleneckVsVector_24Aug2026'
    tmp_xlsx = os.path.join(WORK_DIR, out_name + '_tmp.xlsx')
    out_csv = os.path.join(OUT_DIR, out_name + '.csv')

    log("writing report...")
    smg.Stat_Maker_Gui._write_workbook(None, tmp_xlsx, csv_paths, DATASET_LABELS, None,
                                        bayes_overdispersion, stats, bayesian_stats, bqp_data,
                                        genes, has_bait2, log, ppm_data=ppm_data)
    batch.xlsx_to_csv(tmp_xlsx, out_csv)
    os.remove(tmp_xlsx)
    log("wrote %s" % out_csv)
    log("=== DONE ===")


if __name__ == '__main__':
    main()
