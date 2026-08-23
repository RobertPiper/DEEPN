import os
import re
import csv
import shutil
import pickle
import subprocess


CRAN_DOWNLOAD_URL = 'https://cran.r-project.org/bin/macosx/'
# CRAN's rjags package hard-requires JAGS 4.x (its own configure script
# rejects JAGS 5.x outright: "rjags requires JAGS version 4.x.y"). Homebrew's
# `jags` formula currently only builds 5.0.0, so it is NOT a usable install
# path here - point directly at the official 4.x installers instead.
JAGS_DOWNLOAD_URL = 'https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/'
# Direct download link for the specific installer (matches R 4.3.0+, which
# covers this app's R requirement) - resolves through SourceForge's mirror
# redirect chain to the real .pkg, confirmed working via curl -L.
JAGS_PKG_URL = 'https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/JAGS-4.3.2.pkg/download'
JAGS_PKG_FILENAME = 'JAGS-4.3.2.pkg'

# optparse/dplyr/tidyr/readr/stringr are on CRAN; DESeq2 is Bioconductor-only,
# so installation goes through BiocManager for all of them (it can install
# CRAN packages too).
REQUIRED_R_PACKAGES = ['optparse', 'DESeq2', 'dplyr', 'tidyr', 'readr', 'stringr']

# rjags/runjags are the R-side bridge to the JAGS engine - same tier as
# DESeq2 etc. (installable through R's own package manager), unlike JAGS
# itself, which is a separate standalone program R can't install for you.
# edgeR is deepn's own overdispersion-estimation dependency (its DESCRIPTION
# declares Depends: rjags, runjags, edgeR) - Bioconductor-only, same as DESeq2.
BAYESIAN_R_PACKAGES = ['rjags', 'runjags', 'edgeR']

# The vendored, unmodified copy of pbreheny/deepn (see r_scripts/deepn_bayesian_src/
# README.md for provenance). Installed as a real local-source R package (not
# just sourced loose) so analyzeDeepn()'s own system.file(package="deepn")
# lookups for the JAGS model files resolve correctly, exactly as upstream
# intended - no monkey-patching of the vendored code itself.
DEEPN_PKG_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', 'r_scripts', 'deepn_bayesian_src'))


def check_deepn_package(rscript_exe):
    """Returns True if the vendored 'deepn' R package is installed."""
    proc = subprocess.run(
        [rscript_exe, '-e', 'cat(requireNamespace("deepn", quietly=TRUE))'],
        capture_output=True, text=True
    )
    return proc.stdout.strip() == 'TRUE'


def install_deepn_package(rscript_exe, progress_callback=None):
    """Installs the vendored deepn package from local source (not CRAN/
    Bioconductor - it isn't published there). Requires rjags/runjags/edgeR
    to already be installed (deepn's own Depends). Streams output lines to
    progress_callback(line). Raises RScriptError on failure."""
    expr = 'install.packages(%s, repos=NULL, type="source")' % repr(DEEPN_PKG_DIR).replace("'", '"')
    proc = subprocess.Popen(
        [rscript_exe, '-e', expr],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )
    lines = []
    for line in proc.stdout:
        lines.append(line.rstrip('\n'))
        if progress_callback:
            progress_callback(line.rstrip('\n'))
    proc.wait()
    if proc.returncode != 0:
        raise RScriptError("Installing vendored 'deepn' package failed:\n" + "\n".join(lines[-40:]))
    return "\n".join(lines)


def find_rscript():
    for candidate in ('/usr/local/bin/Rscript', '/opt/homebrew/bin/Rscript'):
        if os.path.exists(candidate):
            return candidate
    found = shutil.which('Rscript')
    if found:
        return found
    return None


def find_jags():
    """Locates a jags executable on PATH/common install locations. Does NOT
    check version compatibility - use check_jags_version() for that, since a
    present-but-wrong-version jags (e.g. Homebrew's 5.x) is not usable by
    rjags and must not be reported as verified."""
    for candidate in ('/usr/local/bin/jags', '/opt/homebrew/bin/jags'):
        if os.path.exists(candidate):
            return candidate
    found = shutil.which('jags')
    if found:
        return found
    return None


def check_jags_version(jags_exe):
    """Returns (compatible, version_string). compatible is True only for
    JAGS 4.x - rjags's own configure script hard-rejects JAGS 5.x, so a 5.x
    install is present but not usable.

    JAGS has no '--version' flag (its CLI treats any argument as a command
    script to open); the version instead appears in the startup banner
    printed to stdout in interactive mode, e.g.
    "Welcome to JAGS 4.3.2 (official binary) on <date>". Feed it an
    immediate 'exit' on stdin so it doesn't hang waiting for commands."""
    try:
        proc = subprocess.run([jags_exe], input='exit\n', capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return False, 'unknown'
    output = (proc.stdout or '') + (proc.stderr or '')
    m = re.search(r'Welcome to JAGS\s+(\d+)\.(\d+)\.(\d+)', output)
    if not m:
        return False, 'unknown'
    version = '.'.join(m.groups())
    return (m.group(1) == '4'), version


def download_jags_pkg(progress_callback=None):
    """Downloads the official JAGS 4.3.2 macOS installer to ~/Downloads via
    curl (follows SourceForge's mirror-redirect chain). Returns the local
    path. Raises RScriptError on failure. Does NOT run/install it - that is
    a separate, explicit step (see open_installer) since it requires the
    user's own admin password in the OS's own installer dialog."""
    dest_dir = os.path.expanduser('~/Downloads')
    dest_path = os.path.join(dest_dir, JAGS_PKG_FILENAME)
    if progress_callback:
        progress_callback("Downloading %s to %s ..." % (JAGS_PKG_FILENAME, dest_dir))
    proc = subprocess.run(
        ['curl', '-L', '--fail', '-sS', '-o', dest_path, JAGS_PKG_URL],
        capture_output=True, text=True
    )
    if proc.returncode != 0 or not os.path.exists(dest_path):
        raise RScriptError("Downloading JAGS installer failed:\n" + (proc.stderr or '').strip())
    if progress_callback:
        size_mb = os.path.getsize(dest_path) / 1e6
        progress_callback("Downloaded %s (%.1f MB)." % (dest_path, size_mb))
    return dest_path


def open_installer(pkg_path):
    """Opens a .pkg in macOS's own Installer.app. The user still has to click
    through it and enter their admin password there - that step cannot be,
    and should not be, automated."""
    subprocess.run(['open', pkg_path])


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


def load_raw_counts_summed(csv_path):
    """Parse a gene_count_summary CSV and reconstruct approximate integer
    read counts from its PPM column: K = round(PPM * TotalReads / 1e6),
    matching the reproducibility-analysis convention (PPM is already
    depth-normalized, so this recovers the pre-normalization count DESeq2
    expects as raw input). Summed across duplicate gene-name rows, same as
    load_ppm_summed(). TotalReads is read from the file's own header line
    (" , TotalReads , <N>", two lines below "File:,...")."""
    total_reads = None
    ppm = {}
    reading = False
    with open(csv_path) as fh:
        for line in fh:
            if total_reads is None and "TotalReads" in line:
                parts = line.split(",")
                if len(parts) >= 3:
                    try:
                        total_reads = float(parts[2])
                    except ValueError:
                        pass
                continue
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
    if total_reads is None:
        raise ValueError("Could not find 'TotalReads' in header of %s" % csv_path)
    return {gene: int(round(v * total_reads / 1e6)) for gene, v in ppm.items()}


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


def build_collated_input_v3(ppm_data, output_dir, collated_path, has_bait2,
                            p_val_thr=0.05, fold_change_thr='none'):
    """v3: same collated-input format as build_collated_input(), but the
    'Bait' group is renamed 'Bait1' (to sit alongside an optional 'Bait2'),
    and a Bait2 CONFIG_START row is only written when has_bait2 is True -
    the R script's group-building loop is already generic over however many
    CONFIG rows exist, so an absent Bait2 row means Bait2 simply doesn't
    exist as a DESeq2 group, and the specificity contrast (Bait1 vs Bait2)
    is skipped rather than failing.

    ppm_data keys: Bait1S, Bait1N, and (if has_bait2) Bait2S, Bait2N, plus
    Vector1S, Vector1N, Vector2S, Vector2N."""
    genes = sorted(ppm_data['Bait1S'].keys())
    with open(collated_path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['data_column', 'note'])
        w.writerow(['ppm', 'DEEPN gene_count_summary PPM export'])
        w.writerow([])
        w.writerow(['CONFIG_START'])
        w.writerow(['enrich_fold_change', 'enrich_p_val', 'normalized', 'output_directory',
                    'sample_names_selected', 'sample_names_background'])
        w.writerow([fold_change_thr, p_val_thr, 'TRUE', output_dir, 'Bait1S', 'Bait1N'])
        if has_bait2:
            w.writerow([fold_change_thr, p_val_thr, 'TRUE', output_dir, 'Bait2S', 'Bait2N'])
        w.writerow([fold_change_thr, p_val_thr, 'TRUE', output_dir, 'Vector1S;Vector2S', 'Vector1N;Vector2N'])
        w.writerow(['DATA_START'])
        data_cols = ['GENE', 'Bait1S', 'Bait1N']
        if has_bait2:
            data_cols += ['Bait2S', 'Bait2N']
        data_cols += ['Vector1S', 'Vector2S', 'Vector1N', 'Vector2N']
        w.writerow(data_cols)
        for g in genes:
            row = [g,
                   round(ppm_data['Bait1S'].get(g, 0.0), 6),
                   round(ppm_data['Bait1N'].get(g, 0.0), 6)]
            if has_bait2:
                row += [round(ppm_data['Bait2S'].get(g, 0.0), 6),
                        round(ppm_data['Bait2N'].get(g, 0.0), 6)]
            row += [round(ppm_data['Vector1S'].get(g, 0.0), 6),
                    round(ppm_data['Vector2S'].get(g, 0.0), 6),
                    round(ppm_data['Vector1N'].get(g, 0.0), 6),
                    round(ppm_data['Vector2N'].get(g, 0.0), 6)]
            w.writerow(row)
    return collated_path


def build_collated_input_v5(raw_count_data, output_dir, collated_path, has_bait2,
                            p_val_thr=0.05, fold_change_thr='none'):
    """v5: same layout as build_collated_input_v3(), but carries raw
    integer counts (from load_raw_counts_summed()) instead of PPM, with
    normalized='FALSE' so run_Y2H_enrichement_stats_v5.R actually computes
    its own size factors (poscounts) instead of treating every sample as
    pre-normalized (size factor 1). PPM was bypassing DESeq2's own
    normalization entirely - this restores it as raw counts, matching the
    only configuration poscounts was validated against.

    raw_count_data keys: Bait1S, Bait1N, and (if has_bait2) Bait2S, Bait2N,
    plus Vector1S, Vector1N, Vector2S, Vector2N - same shape as
    build_collated_input_v3's ppm_data, just integer counts instead of PPM."""
    genes = sorted(raw_count_data['Bait1S'].keys())
    with open(collated_path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['data_column', 'note'])
        w.writerow(['raw_count', 'Reconstructed from DEEPN gene_count_summary PPM (K = round(PPM * TotalReads / 1e6))'])
        w.writerow([])
        w.writerow(['CONFIG_START'])
        w.writerow(['enrich_fold_change', 'enrich_p_val', 'normalized', 'output_directory',
                    'sample_names_selected', 'sample_names_background'])
        w.writerow([fold_change_thr, p_val_thr, 'FALSE', output_dir, 'Bait1S', 'Bait1N'])
        if has_bait2:
            w.writerow([fold_change_thr, p_val_thr, 'FALSE', output_dir, 'Bait2S', 'Bait2N'])
        w.writerow([fold_change_thr, p_val_thr, 'FALSE', output_dir, 'Vector1S;Vector2S', 'Vector1N;Vector2N'])
        w.writerow(['DATA_START'])
        data_cols = ['GENE', 'Bait1S', 'Bait1N']
        if has_bait2:
            data_cols += ['Bait2S', 'Bait2N']
        data_cols += ['Vector1S', 'Vector2S', 'Vector1N', 'Vector2N']
        w.writerow(data_cols)
        for g in genes:
            row = [g,
                   raw_count_data['Bait1S'].get(g, 0),
                   raw_count_data['Bait1N'].get(g, 0)]
            if has_bait2:
                row += [raw_count_data['Bait2S'].get(g, 0),
                        raw_count_data['Bait2N'].get(g, 0)]
            row += [raw_count_data['Vector1S'].get(g, 0),
                    raw_count_data['Vector2S'].get(g, 0),
                    raw_count_data['Vector1N'].get(g, 0),
                    raw_count_data['Vector2N'].get(g, 0)]
            w.writerow(row)
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


def run_bayesian_r_script(r_script_path, csv_paths, threshold, out_dir, rscript_exe=None):
    """Runs the Bayesian/MCMC (analyzeDeepn/JAGS) R script directly against the
    raw gene_count_summary.csv paths already collected by the GUI (Vector1S/
    Vector1N/Vector2S/Vector2N/Bait1S/Bait1N and optionally Bait2S/Bait2N) -
    the vendored deepn package reads this same file format itself, so no
    intermediate collation step is needed here (unlike the DESeq2/PPM side).
    Returns (stdout_text, outfile_path)."""
    exe = rscript_exe or find_rscript()
    if exe is None:
        raise RScriptError("Could not find Rscript on this machine. Install R (e.g. via Homebrew: 'brew install r').")

    outfile = os.path.join(out_dir, 'bayesian_stats.csv')
    msgfile = os.path.join(out_dir, 'bayesian_overdispersion.txt')
    args = [
        exe, r_script_path,
        '--vec_sel1', csv_paths['Vector1S'], '--vec_nonsel1', csv_paths['Vector1N'],
        '--vec_sel2', csv_paths['Vector2S'], '--vec_nonsel2', csv_paths['Vector2N'],
        '--bait1_sel', csv_paths['Bait1S'], '--bait1_nonsel', csv_paths['Bait1N'],
        '--threshold', str(threshold),
        '--outfile', outfile, '--msgfile', msgfile,
    ]
    if 'Bait2S' in csv_paths and 'Bait2N' in csv_paths:
        args += ['--bait2_sel', csv_paths['Bait2S'], '--bait2_nonsel', csv_paths['Bait2N']]

    proc = subprocess.run(args, capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RScriptError("Bayesian R script failed (exit %d):\n%s" % (proc.returncode, output))

    # analyzeDeepn()'s own overdispersion report (edgeR common dispersion,
    # computed on raw counts - the original Bayesian-method math, distinct
    # from DESeq2's own dispersion estimate on the PPM side).
    bayesian_overdispersion = {}
    patterns = {
        'vector_nonselect': r'Baseline \(vector only\):\s*([0-9.eE+-]+)',
        'vector_bait_nonselect': r'Baseline \(vector \+ bait\):\s*([0-9.eE+-]+)',
        'vector_select': r'Selection:\s*([0-9.eE+-]+)',
    }
    for key, pattern in patterns.items():
        m = re.search(pattern, output)
        if m:
            try:
                bayesian_overdispersion[key] = float(m.group(1))
            except ValueError:
                pass

    return output, outfile, bayesian_overdispersion


def find_bqp_path(directory, basename):
    """Locates the 5p .bqp file for a gene_count_summary basename, checking
    both naming conventions: the current 5p_<basename>.bqp, and the legacy
    convention (pre-dating the 5p_/3p_ prefix scheme) where 5p files were
    written with no prefix at all and only 3p files got a '3p_' prefix - a
    plain <basename>.bqp sitting next to a 3p_<basename>.bqp sibling is the
    legacy 5p file. Confirmed against real 2018-era data (DEEPN_2018_RabGTPase),
    which has zero 5p_-prefixed files but a full set of unprefixed ones.
    Returns the resolved path, or None if neither exists."""
    query_dir = os.path.join(directory, 'blast_results_query')
    prefixed = os.path.join(query_dir, '5p_' + basename + '.bqp')
    if os.path.exists(prefixed):
        return prefixed
    legacy = os.path.join(query_dir, basename + '.bqp')
    if os.path.exists(legacy):
        return legacy
    return None


def load_bayesian_stats_output(bayesian_stats_csv):
    """Reads the Bayesian R script's output CSV into dict[gene] = row_dict."""
    stats = {}
    with open(bayesian_stats_csv, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gene = row.get('Gene')
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
            # round(), not format(...,'.1f') - format() returns a string,
            # which Excel then treats as text rather than a sortable number.
            gene_stat[key] = round(gene_stat[key] * 100.0 / gene_stat['total'], 1)
    except ZeroDivisionError:
        for key in ('frame_orf', 'upstream', 'in_orf', 'downstream', 'in_frame', 'backwards', 'intron'):
            gene_stat[key] = 0
    return gene_stat
