#!/usr/bin/env python3
"""Rebuilds Rab_unified_criteria.xlsx with the final 62-gene garbage/
background union list now subtracted from the new-hit calls, distinguishing
genes lost specifically because they're now recognized as garbage from
genes that simply never cleared the hit criteria at all.
"""
import sys, os, re, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/Users/robertpiper2/LOCAL_DEEPN')
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
os.chdir('/Users/robertpiper2/LOCAL_DEEPN')

import xlsxwriter as xls
import stat_maker_gui_v7 as smg
import sm6_batch as batch
import rab_hit_recapitulation as rhr

OUT_PATH = '/Users/robertpiper2/endosomeLAB Dropbox/Rob Piper/SM6_batch_output/Rab_unified_criteria.xlsx'
GARBAGE_CSV = '/Users/robertpiper2/Dropbox (endosomeLAB)/SM6_batch_output/Rab_garbage_union_final.csv'


def construct_label(filename, expect):
    m = re.search(r'(QL|TN|SN)', filename)
    return m.group(1) if m else expect


def main():
    print("Loading GUI defaults, reference hits, SM6 data, and garbage union list...")
    win = smg.Stat_Maker_Gui()
    CRITERIA = win._criteria_dict(win.b1_criteria)
    ref, col_meta = rhr.load_reference_hits()
    data = rhr.load_all_sm6_data()
    GARBAGE_UNION = {r['Gene'] for r in csv.DictReader(open(GARBAGE_CSV))}
    print("Garbage union list:", len(GARBAGE_UNION), "genes")

    construct = {}
    for fn in batch.CANONICAL:
        rid = batch.rab_id(fn)
        cfg = batch.get_config(os.path.join(batch.SM, fn))
        construct[(rid, 1)] = construct_label(cfg['Bait1_Selected_1'], 'QL')
        construct[(rid, 2)] = construct_label(cfg['Bait2_Selected_1'], 'TN')

    print("Computing new hit sets, with garbage subtracted...")
    col_data = {}
    all_new_genes = set()
    total_gained_inst = total_lost_inst = total_retained_inst = 0
    total_lost_to_garbage_inst = 0
    for _c, rid, bait_num in col_meta:
        maps = data[rid][bait_num]
        genes = data[rid]['genes']
        ref_set = ref.get(rid, {}).get(bait_num, set())
        passes_criteria = {g for g in genes if smg._gene_passes_criteria(g, CRITERIA, maps)}
        new_hits = passes_criteria - GARBAGE_UNION
        gained = new_hits - ref_set
        lost = ref_set - new_hits
        retained = new_hits & ref_set
        lost_to_garbage = lost & GARBAGE_UNION
        lost_to_criteria = lost - lost_to_garbage
        col_data[(rid, bait_num)] = {'new_hits': new_hits, 'gained': gained, 'lost': lost,
                                      'retained': retained, 'lost_to_garbage': lost_to_garbage,
                                      'lost_to_criteria': lost_to_criteria}
        all_new_genes |= new_hits
        total_gained_inst += len(gained); total_lost_inst += len(lost); total_retained_inst += len(retained)
        total_lost_to_garbage_inst += len(lost_to_garbage)

    total_ref = sum(len(ref.get(rid, {}).get(bait_num, set())) for _c, rid, bait_num in col_meta)
    distinct_gained = set().union(*(v['gained'] for v in col_data.values()))
    distinct_lost = set().union(*(v['lost'] for v in col_data.values()))
    distinct_lost_to_garbage = set().union(*(v['lost_to_garbage'] for v in col_data.values()))
    print("total prior hits: %d | retained: %d (%.1f%%) | lost: %d instances / %d genes "
          "(of which %d instances / %d genes lost specifically to garbage) | gained: %d instances / %d genes"
          % (total_ref, total_retained_inst, 100.0*total_retained_inst/total_ref, total_lost_inst, len(distinct_lost),
             total_lost_to_garbage_inst, len(distinct_lost_to_garbage), total_gained_inst, len(distinct_gained)))

    print("Writing workbook...")
    wb = xls.Workbook(OUT_PATH)
    ws = wb.add_worksheet('Unified Criteria')

    fmt_title = wb.add_format({'bold': True, 'font_size': 14})
    fmt_bold = wb.add_format({'bold': True})
    fmt_header = wb.add_format({'bold': True, 'text_wrap': True, 'align': 'center', 'valign': 'vcenter',
                                 'bg_color': '#D9E8FB', 'border': 1})
    fmt_gene = wb.add_format({'border': 1})
    fmt_gained = wb.add_format({'border': 1, 'bg_color': '#C6EFCE'})        # green
    fmt_lost = wb.add_format({'border': 1, 'bg_color': '#FFC7CE'})         # pink - lost, criteria not met
    fmt_lost_garbage = wb.add_format({'border': 1, 'bg_color': '#FFD966'})  # gold - lost, now flagged garbage
    fmt_hit_green = wb.add_format({'border': 1, 'bg_color': '#C6EFCE', 'align': 'center'})
    fmt_section = wb.add_format({'bold': True, 'italic': True})

    row = 0
    ws.write(row, 0, 'Rab_unified_criteria — Revised Rab GTPase Interactor Table (Stat Maker v7)', fmt_title)
    row += 2
    ws.write(row, 0, 'Default hit criteria (locked in GUI, 2026-08-24):', fmt_bold)
    row += 1
    ws.write(row, 0, 'in-frame:Forward >= 80%   |   AdjEnr (Bait vs Vector) > 1x   |   '
                      'DESeq2 log2FC (Bait vs Vector) > 3x   |   p-value_raw <= 0.01   |   '
                      'P (MCMC) > 0.8   |   p-value_norm: not applied (unreliable outside low-crush datasets)')
    row += 1
    ws.write(row, 0, 'Garbage/background union list applied (%d genes: old Word-doc list + vector-only DESeq2 '
                      '>=8/16 + bottleneck DESeq2 >=2/4 reference comparisons):' % len(GARBAGE_UNION), fmt_bold)
    row += 2
    ws.write(row, 0, 'Compared to the prior manually-curated hit list (819 hits, 72 Rab-bait datasets):', fmt_bold)
    row += 1
    ws.write(row, 0, 'Retained: %d (%.1f%%)   |   Lost: %d instances / %d distinct genes   |   '
                      'Gained: %d instances / %d distinct genes (shaded green below)'
             % (total_retained_inst, 100.0*total_retained_inst/total_ref, total_lost_inst, len(distinct_lost),
                total_gained_inst, len(distinct_gained)))
    row += 1
    ws.write(row, 0, 'Of the lost genes: %d instances / %d distinct genes lost specifically because they are now '
                      'on the garbage list (shaded gold below); the rest simply did not clear the hit criteria '
                      '(shaded pink below).' % (total_lost_to_garbage_inst, len(distinct_lost_to_garbage)))
    row += 2
    ws.write(row, 0, 'LEFT: per-Rab-bait verified hit lists under the new criteria, garbage-subtracted '
                      '(green = newly gained; below each column, prior hits NOT recaptured - gold if now flagged '
                      'garbage, pink if the criteria simply weren\'t met).   '
                      'RIGHT: master gene list x Rab-bait matrix, green HIT cells (mirrors the original table layout).',
             fmt_section)
    row += 2

    header_row = row
    data_start = header_row + 1

    max_new = max(len(v['new_hits']) for v in col_data.values())
    max_lost = max(len(v['lost']) for v in col_data.values())
    for col, (_c, rid, bait_num) in enumerate(col_meta):
        label = '%s_B%d (%s)' % (rid, bait_num, construct[(rid, bait_num)])
        ws.write(header_row, col, label, fmt_header)
        ws.set_column(col, col, 13)
        d = col_data[(rid, bait_num)]
        new_sorted = sorted(d['new_hits'])
        for i, g in enumerate(new_sorted):
            fmt = fmt_gained if g in d['gained'] else fmt_gene
            ws.write(data_start + i, col, g, fmt)

        lost_label_row = data_start + max_new + 1
        lost_start = lost_label_row + 1
        if col == 0:
            ws.write(lost_label_row, 0, 'LOST (prior hit, not recaptured - gold=now garbage, pink=criteria not met):',
                     fmt_section)
        lost_sorted = sorted(d['lost'])
        for i, g in enumerate(lost_sorted):
            fmt = fmt_lost_garbage if g in d['lost_to_garbage'] else fmt_lost
            ws.write(lost_start + i, col, g, fmt)

    left_block_end_col = len(col_meta) - 1

    right_start_col = left_block_end_col + 3
    gene_col = right_start_col
    ws.write(header_row, gene_col, 'Gene', fmt_header)
    ws.set_column(gene_col, gene_col, 14)
    for i, (_c, rid, bait_num) in enumerate(col_meta):
        col = right_start_col + 1 + i
        label = '%s_B%d (%s)' % (rid, bait_num, construct[(rid, bait_num)])
        ws.write(header_row, col, label, fmt_header)
        ws.set_column(col, col, 13)

    for i, g in enumerate(sorted(all_new_genes)):
        r = data_start + i
        ws.write(r, gene_col, g, fmt_bold)
        for j, (_c, rid, bait_num) in enumerate(col_meta):
            col = right_start_col + 1 + j
            if g in col_data[(rid, bait_num)]['new_hits']:
                is_gained = g in col_data[(rid, bait_num)]['gained']
                ws.write(r, col, 'HIT', fmt_hit_green if is_gained else fmt_gene)

    ws.freeze_panes(data_start, 1)
    wb.close()
    print("Wrote %s" % OUT_PATH)
    print("LEFT block: %d columns, up to %d new hits + %d lost per column" % (len(col_meta), max_new, max_lost))
    print("RIGHT block: %d gene rows x %d Rab-bait columns" % (len(all_new_genes), len(col_meta)))


if __name__ == '__main__':
    main()
