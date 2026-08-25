#!/usr/bin/env python3
"""Runs the bottleneck-pool-vs-vector Stat Maker v7 pipeline 4 more times,
each against a different Vector1/Vector2 reference pair, to test how
sensitive the bottleneck hit list is to which specific vector prep is used
as reference. For each run, the bottleneck's own Bait1-Non-Selected partner
is a DIFFERENT vector, never one of that run's own Vector1/Vector2
reference pair, keeping all data independent (this is why the original
Vector5+Vector6 request was adjusted - Vector5 was both the reference and
the bottleneck's own Bait1N).
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
WORK_ROOT = '/private/tmp/claude-503/-Users-robertpiper2-Desktop/d0bd16a2-ab29-49c2-9981-9811e0c1149f/scratchpad/bottleneck_multi_ref'
LOG_PATH = os.path.join(WORK_ROOT, 'log.txt')
BAIT1S = 'HuORF_pTEF_Control_S202_SEL_summary.csv'

RUNS = [
    {'name': 'Vector13_Vector4', 'bait1n': 'pTEFGBD_Vector3_S19_NON_summary.csv',
     'v1s': 'pTEFGBD_Vector13_S125_SEL_summary.csv', 'v1n': 'pTEFGBD_Vector13_S124_NON_summary.csv',
     'v2s': 'pTEFGBD_Vector4_S38_SEL_summary.csv', 'v2n': 'pTEFGBD_Vector4_S37_NON_summary.csv'},
    {'name': 'Vector5_Vector6', 'bait1n': 'pTEFGBD_Vector7_S67_NON_summary.csv',
     'v1s': 'pTEFGBD_Vector5_S48_SEL_summary.csv', 'v1n': 'pTEFGBD_Vector5_S47_NON_summary.csv',
     'v2s': 'pTEFGBD_Vector6_S58_SEL_summary.csv', 'v2n': 'pTEFGBD_Vector6_S57_NON_summary.csv'},
    {'name': 'Vector9_Vector1plain', 'bait1n': 'pTEFGBD_Vector10_S139_NON_summary.csv',
     'v1s': 'pTEFGBD_Vector9_S96_SEL_summary.csv', 'v1n': 'pTEFGBD_Vector9_S95_NON_summary.csv',
     'v2s': 'pTEFGBDvector1_S8_SEL_summary.csv', 'v2n': 'pTEFGBDvector1_S7_NON_summary.csv'},
    {'name': 'Vector12_Vector8', 'bait1n': 'pTEFstarGBD_Vector1_S182_NON_summary.csv',
     'v1s': 'pTEFGBD_Vector12_S111_SEL_summary.csv', 'v1n': 'pTEFGBD_Vector12_S110_NON_summary.csv',
     'v2s': 'pTEFGBD_Vector8_S82_SEL_summary.csv', 'v2n': 'pTEFGBD_Vector8_S81_NON_summary.csv'},
]


def log(msg):
    line = "[%s] %s" % (time.strftime('%H:%M:%S'), msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a') as f:
        f.write(line + '\n')


def run_one(cfg, rscript_exe):
    name = cfg['name']
    work_dir = os.path.join(WORK_ROOT, name)
    os.makedirs(work_dir, exist_ok=True)
    csv_paths = {
        'Bait1S': GC + BAIT1S, 'Bait1N': GC + cfg['bait1n'],
        'Vector1S': GC + cfg['v1s'], 'Vector1N': GC + cfg['v1n'],
        'Vector2S': GC + cfg['v2s'], 'Vector2N': GC + cfg['v2n'],
    }
    for k, p in csv_paths.items():
        if not os.path.exists(p):
            raise FileNotFoundError("%s: %s" % (k, p))
    has_bait2 = False
    DATASET_LABELS = {
        'Bait1N': 'Bait1_Non-Selected (%s)' % cfg['bait1n'], 'Bait1S': 'Bait1_Selected (Bottleneck pool)',
        'Vector1N': 'Vector_Non-Selected_1', 'Vector2N': 'Vector_Non-Selected_2',
        'Vector1S': 'Vector_Selected_1', 'Vector2S': 'Vector_Selected_2',
    }

    log("[%s] loading raw counts..." % name)
    raw_count_data = {k: sc.load_raw_counts_summed(p) for k, p in csv_paths.items()}
    genes = sorted(raw_count_data['Bait1S'].keys())
    ppm_data = {k: sc.load_ppm_summed(p) for k, p in csv_paths.items()}

    collated_path = os.path.join(work_dir, 'collated_input.csv')
    sc.build_collated_input_v5(raw_count_data, work_dir, collated_path, has_bait2)
    log("[%s] running raw-counts + poscounts DESeq2..." % name)
    sc.run_r_script(smg.R_SCRIPT_PATH, collated_path, rscript_exe)
    stats = sc.load_stats_output(os.path.join(work_dir, 'everything_combined.csv'))
    stats = smg._suffix_bait_vs_vector(stats, has_bait2, '_raw')

    log("[%s] running normalized-PPM DESeq2..." % name)
    norm_dir = os.path.join(work_dir, 'norm')
    os.makedirs(norm_dir, exist_ok=True)
    collated_norm = os.path.join(norm_dir, 'collated_input.csv')
    sc.build_collated_input_v3(ppm_data, norm_dir, collated_norm, has_bait2)
    sc.run_r_script(smg.R_SCRIPT_PATH, collated_norm, rscript_exe)
    stats_norm = sc.load_stats_output(os.path.join(norm_dir, 'everything_combined.csv'))
    smg._merge_bait_vs_vector_norm(stats, stats_norm, has_bait2)

    log("[%s] running Bayesian/MCMC..." % name)
    _, bayes_csv, bayes_overdispersion = sc.run_bayesian_r_script(
        smg.BAYESIAN_R_SCRIPT_PATH, csv_paths, smg.BAYESIAN_THRESHOLD_PPM, work_dir, rscript_exe)
    bayesian_stats = sc.load_bayesian_stats_output(bayes_csv)

    log("[%s] loading junction data..." % name)
    bqp_data = {}
    for k, p in csv_paths.items():
        basename = os.path.basename(p).replace('_summary.csv', '')
        bqp_path = sc.find_bqp_path(BQP_ROOT, basename)
        bqp_data[k] = sc.load_bqp(bqp_path) if bqp_path else {}

    out_csv = os.path.join(OUT_DIR, 'SM6_Bottleneck_%s_24Aug2026.csv' % name)
    tmp_xlsx = os.path.join(work_dir, 'tmp.xlsx')
    log("[%s] writing report..." % name)
    smg.Stat_Maker_Gui._write_workbook(None, tmp_xlsx, csv_paths, DATASET_LABELS, None,
                                        bayes_overdispersion, stats, bayesian_stats, bqp_data,
                                        genes, has_bait2, log, ppm_data=ppm_data)
    batch.xlsx_to_csv(tmp_xlsx, out_csv)
    os.remove(tmp_xlsx)
    log("[%s] DONE -> %s" % (name, out_csv))
    return out_csv


def main():
    os.makedirs(WORK_ROOT, exist_ok=True)
    open(LOG_PATH, 'w').close()
    rscript_exe = sc.find_rscript()
    for cfg in RUNS:
        run_one(cfg, rscript_exe)
    log("=== ALL 4 RUNS DONE ===")


if __name__ == '__main__':
    main()
