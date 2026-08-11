import os
import re
import csv
import shutil
import pickle
import subprocess


CRAN_DOWNLOAD_URL = 'https://cran.r-project.org/bin/macosx/'

# optparse/dplyr/tidyr/readr/stringr are on CRAN; DESeq2 is Bioconductor-only,
# so installation goes through BiocManager for all of them (it can install
# CRAN packages too).
REQUIRED_R_PACKAGES = ['optparse', 'DESeq2', 'dplyr', 'tidyr', 'readr', 'stringr']


def find_rscript():
    for candidate in ('/usr/local/bin/Rscript', '/opt/homebrew/bin/Rscript'):
        if os.path.exists(candidate):
            return candidate
    found = shutil.which('Rscript')
    if found:
        return found
    return None


def check_r_packages(rscript_exe, packages=None):
    """Returns dict[package_name] = True/False (installed or not)."""
    packages = packages or REQUIRED_R_PACKAGES
    check_expr = ";".join(
        'cat("%s:", requireNamespace("%s", quietly=TRUE), "\\n")' % (p, p)
        for p in packages
    )
    proc = subprocess.run([rscript_exe, '-e', check_expr], capture_output=True, text=True)
    status = {}
    for line in proc.stdout.splitlines():
        m = re.match(r'^([A-Za-z0-9._]+):\s*(TRUE|FALSE)\s*$', line.strip())
        if m:
            status[m.group(1)] = (m.group(2) == 'TRUE')
    for p in packages:
        status.setdefault(p, False)
    return status


def install_r_packages(rscript_exe, packages, progress_callback=None):
    """Installs R packages via BiocManager (works for both CRAN and Bioconductor
    packages). Streams output lines to progress_callback(line) as they arrive.
    Raises RScriptError on failure."""
    pkg_vector = "c(" + ", ".join('"%s"' % p for p in packages) + ")"
    install_expr = (
        'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); '
        'BiocManager::install(%s, update=FALSE, ask=FALSE)' % pkg_vector
    )
    proc = subprocess.Popen(
        [rscript_exe, '-e', install_expr],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )
    lines = []
    for line in proc.stdout:
        lines.append(line.rstrip('\n'))
        if progress_callback:
            progress_callback(line.rstrip('\n'))
    proc.wait()
    if proc.returncode != 0:
        raise RScriptError("R package installation failed:\n" + "\n".join(lines[-40:]))
    return "\n".join(lines)


def load_ppm_summed(csv_path):
    """Parse a gene_count_summary CSV, summing PPM across any duplicate
    gene-name rows (genes annotated on multiple alt/scaffold contigs)."""
    ppm = {}
    reading = False
    with open(csv_path) as fh:
        for line in fh:
            if line.startswith("Chromosome"):
                reading = True
                continue
            if not reading:
                continue
            parts = line.split(",")
            if len(parts) < 3:
                continue
            gene = parts[1].strip()
            try:
                val = float(parts[2])
            except ValueError:
                continue
            ppm[gene] = ppm.get(gene, 0.0) + val
    return ppm


def count_genes(csv_path):
    return len(load_ppm_summed(csv_path))


def check_gene_counts_consistent(csv_paths):
    """Returns (consistent: bool, counts: dict[path -> count])."""
    counts = {p: count_genes(p) for p in csv_paths}
    consistent = len(set(counts.values())) <= 1
    return consistent, counts


def build_collated_input(ppm_data, output_dir, collated_path, p_val_thr=0.05, fold_change_thr='none'):
    """ppm_data: dict with keys BaitS, BaitN, Vector1S, Vector1N, Vector2S, Vector2N
    -> dict[gene_name] = ppm (float)."""
    genes = sorted(ppm_data['BaitS'].keys())
    with open(collated_path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['data_column', 'note'])
        w.writerow(['ppm', 'DEEPN gene_count_summary PPM export'])
        w.writerow([])
        w.writerow(['CONFIG_START'])
        w.writerow(['enrich_fold_change', 'enrich_p_val', 'normalized', 'output_directory',
                    'sample_names_selected', 'sample_names_background'])
        w.writerow([fold_change_thr, p_val_thr, 'TRUE', output_dir, 'BaitS', 'BaitN'])
        w.writerow([fold_change_thr, p_val_thr, 'TRUE', output_dir, 'Vector1S;Vector2S', 'Vector1N;Vector2N'])
        w.writerow(['DATA_START'])
        w.writerow(['GENE', 'BaitS', 'BaitN', 'Vector1S', 'Vector2S', 'Vector1N', 'Vector2N'])
        for g in genes:
            w.writerow([
                g,
                round(ppm_data['BaitS'].get(g, 0.0), 6),
                round(ppm_data['BaitN'].get(g, 0.0), 6),
                round(ppm_data['Vector1S'].get(g, 0.0), 6),
                round(ppm_data['Vector2S'].get(g, 0.0), 6),
                round(ppm_data['Vector1N'].get(g, 0.0), 6),
                round(ppm_data['Vector2N'].get(g, 0.0), 6),
            ])
    return collated_path


class RScriptError(Exception):
    pass


def run_r_script(r_script_path, collated_path, rscript_exe=None):
    """Runs the DESeq2 R script. Returns (stdout_text, nonselect_overdispersion or None)."""
    exe = rscript_exe or find_rscript()
    if exe is None:
        raise RScriptError("Could not find Rscript on this machine. Install R (e.g. via Homebrew: 'brew install r').")

    proc = subprocess.run(
        [exe, r_script_path, '--collated_file', collated_path],
        capture_output=True, text=True
    )
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RScriptError("R script failed (exit %d):\n%s" % (proc.returncode, output))

    overdispersion = None
    match = re.search(r'non-selected samples:\s*([0-9.eE+-]+)', output)
    if match:
        try:
            overdispersion = float(match.group(1))
        except ValueError:
            overdispersion = None
    return output, overdispersion


def load_stats_output(everything_combined_csv):
    """Reads the R script's everything_combined.csv into dict[gene] = row_dict."""
    stats = {}
    with open(everything_combined_csv, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gene = row.get('gene')
            if gene:
                stats[gene] = row
    return stats


def load_bqp(path):
    with open(path, 'rb') as fh:
        return pickle.load(fh)


def get_junction_stats(bqp_data, gene_name):
    """Ported from stat_maker_gui.py's Stat_Maker_Gui.get_junction_stats - computes
    junction-frame/ORF percentages for one gene from a loaded .bqp dict."""
    gene_stat = {'frame_orf': 0, 'upstream': 0, 'in_orf': 0, 'downstream': 0,
                 'in_frame': 0, 'backwards': 0, 'intron': 0, 'total': 0}
    try:
        gene = bqp_data[gene_name]
        for nm in gene.keys():
            if nm == 'stats':
                continue
            for j in gene[nm]:
                if j.orf == 'in_orf' and j.frame == 'in_frame':
                    gene_stat['frame_orf'] += j.ppm
                if j.orf == 'in_orf':
                    gene_stat['in_orf'] += j.ppm
                elif j.orf == 'upstream':
                    gene_stat['upstream'] += j.ppm
                elif j.orf == 'downstream':
                    gene_stat['downstream'] += j.ppm
                if j.frame == 'in_frame':
                    gene_stat['in_frame'] += j.ppm
                elif j.frame == 'backwards':
                    gene_stat['backwards'] += j.ppm
                elif j.frame == 'intron':
                    gene_stat['intron'] += j.ppm
                gene_stat['total'] += j.ppm
    except KeyError:
        pass
    try:
        for key in ('frame_orf', 'upstream', 'in_orf', 'downstream', 'in_frame', 'backwards', 'intron'):
            gene_stat[key] = format(gene_stat[key] * 100.0 / gene_stat['total'], ".1f")
    except ZeroDivisionError:
        for key in ('frame_orf', 'upstream', 'in_orf', 'downstream', 'in_frame', 'backwards', 'intron'):
            gene_stat[key] = '0'
    return gene_stat
