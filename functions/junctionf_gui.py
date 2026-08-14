import os
import sys
import time
import struct
import platform
import pickle
import subprocess
from sys import platform as _platform
from collections import Counter

import threading
import multiprocessing
import joblib
import functions.process as process
import functions.structures as sts

def threaded(fn):
    def wrapper(*args, **kwargs):
        thread = threading.Thread(target=fn, args=args, kwargs=kwargs)
        thread.start()
        return thread
    return wrapper

class junctionf():
    def __init__(self, f, p):
        self.fileio = f
        self.printio = p
        self.process = process.process()
        self.blast_pipe = None

    def sigterm_handler(self, _signo, _stack_frame):
        # blastn is a subprocess of this process, not of DEEPN's launcher -
        # SIGTERM sent to us (via the launcher's Abort button, or DEEPN
        # quitting while this is running) doesn't automatically reach our
        # own child. Without this handler blastn just keeps running,
        # orphaned, using full CPU, with no way to stop it short of killing
        # it manually in Activity Monitor.
        if self.blast_pipe and self.blast_pipe.poll() is None:
            self.blast_pipe.terminate()
            print(">>> Terminated BLAST process. Exiting.")
        else:
            print(">>> Exiting.")
        sys.stdout.flush()
        sys.exit()

    def _getjunction(self, begin):
        var = ''
        infile = open(os.path.join(os.curdir, "Y2Hreadme.txt"), 'r')
        for line in infile:
            if begin in line:
                var = (next(infile))
        infile.close()
        return var

    def _get_unprocessed_files(self, list, suffix1, processed_list, suffix2):
        for processed_file in processed_list:
            for f in list:
                if f.replace(suffix1, "") == processed_file.replace(suffix2, ""):
                    list.remove(f)
                    break
        return list

    _DIRECTION_LABELS = {'5p_': "5'", '3p_': "3'"}

    def junction_search(self, directory, junction_folder, input_data_folder, blast_results_folder,
                        blast_results_query, junctions_array, exclusion_sequence, prefix):
        exclusion_sequence.upper()
        jseqs = self._make_search_junctions(junctions_array)
        direction_label = self._DIRECTION_LABELS.get(prefix, prefix)
        print(">>> Searching for %s junctions. Primary, secondary, and tertiary sequences searched:" % direction_label)
        for seq in jseqs:
            print("      %s" % seq)
        sys.stdout.flush()
        unmapped_sam_files = self.fileio.get_sam_filelist(directory, input_data_folder)

        print('>>> Starting junction search.')
        sys.stdout.flush()

        for f in unmapped_sam_files:
            print('>>> [%s junctions] File: %s' % (direction_label, f))
            sys.stdout.flush()
            filename = os.path.join(directory, input_data_folder, f)
            input_filehandle = open(filename)
            input_file_size = os.path.getsize(filename)
            output_filehandle = open(os.path.join(directory, junction_folder, prefix + f.replace(".sam", '.junctions.txt')), 'w')
            self._search_for_junctions(input_filehandle, jseqs, exclusion_sequence, output_filehandle, f, input_file_size)
            output_filehandle.close()
            input_filehandle.close()
        self._multi_convert(directory, junction_folder, blast_results_folder, prefix)

    def _multi_convert(self, directory, infolder, outfolder, prefix):
        file_list = [f for f in self.fileio.get_file_list(directory, infolder, ".txt") if f.startswith(prefix)]

        print(' ')
        for f in file_list:
            self.fileio.make_FASTA(os.path.join(directory, infolder, f),
                                   os.path.join(directory, outfolder, f[:-4] + ".fa"))

    def _count_fasta_sequences(self, fa_path):
        count = 0
        with open(fa_path, 'r') as fh:
            for line in fh:
                if line.startswith('>'):
                    count += 1
        return count

    def _monitor_blast_progress(self, output_file, total_sequences, stop_event, checkpoint=10000):
        """blastn has no progress flag, but with -outfmt 7 it writes a
        "# Query: ..." line to the output file the moment each query
        finishes - so the growing output file itself is the progress
        signal. Polls every couple seconds (pure I/O wait, negligible CPU)
        and prints an update every `checkpoint` completed sequences, with
        an ETA based on the actual rate over the last checkpoint interval."""
        last_pos = 0
        completed = 0
        last_checkpoint_count = 0
        last_checkpoint_time = time.time()
        while not stop_event.is_set():
            if os.path.exists(output_file):
                with open(output_file, 'r') as fh:
                    fh.seek(last_pos)
                    chunk = fh.read()
                    last_pos = fh.tell()
                completed += chunk.count('# Query:')
                if completed - last_checkpoint_count >= checkpoint:
                    now = time.time()
                    elapsed = now - last_checkpoint_time
                    done_this_round = completed - last_checkpoint_count
                    rate = done_this_round / elapsed if elapsed > 0 else 0
                    remaining = max(0, total_sequences - completed)
                    eta_minutes = (remaining / rate / 60.0) if rate > 0 else 0
                    pct = (completed * 100.0 / total_sequences) if total_sequences else 0
                    # \r, not a trailing newline, same convention as the
                    # junction-scan progress above - updates the current
                    # line in place (both in a real terminal and in deepn.py's
                    # status box, which strips the \r and overwrites the
                    # current line for anything that isn't a new section
                    # header - see start_match in deepn.py) instead of
                    # stacking a new line per checkpoint.
                    sys.stdout.write('\rBLAST progress: %d / %d sequences (%.0f%%) - last %d took %.0fs, ~%.1f min remaining' %
                         (completed, total_sequences, pct, done_this_round, elapsed, eta_minutes))
                    sys.stdout.flush()
                    last_checkpoint_count = completed
                    last_checkpoint_time = now
            stop_event.wait(2.0)

    def blast_search(self, directory, db_name, blast_results_folder, blast_results_query, prefix):
        platform_specific_path = 'osx'
        suffix = ''
        if _platform == "linux" or _platform == "linux2":
            platform_specific_path = 'linux'
        elif _platform == "darwin":
            platform_specific_path = 'osx'
        elif _platform.startswith('win'):
            platform_specific_path = 'windows'
            suffix = '.exe'
        machine = platform.machine()
        if machine == 'arm64':
            arch_dir = 'arm64'
        else:
            arch_dir = 'x' + str(struct.calcsize("P") * 8)
        blast_path = os.path.join(os.curdir, 'ncbi_blast', 'bin', platform_specific_path, arch_dir)
        blast_db = os.path.join(os.curdir, 'ncbi_blast', 'db')
        db_path = os.path.join(blast_db, db_name)
        print(">>> Selected Blast DB: %s" % db_name)
        sys.stdout.flush()
        file_list = [f for f in self.fileio.get_file_list(directory, blast_results_folder, ".fa") if f.startswith(prefix)]

        # Leave one core free rather than requesting every core - otherwise
        # blastn saturates the whole machine and everything else (including
        # this app's own GUI) becomes unresponsive for the run's duration.
        # On macOS there's no API for a regular app to pin threads to
        # specific cores, or even to just the Performance cores - only a
        # thread *count*, which is all this controls; which physical cores
        # actually run those threads is entirely up to the OS scheduler.
        # Windows and Linux both expose explicit CPU affinity to user
        # processes and can report which cores are P-cores vs E-cores, so
        # if this is ever ported there, it may be worth explicitly keeping
        # blastn off the E-cores instead of just reducing the thread count -
        # but that needs testing against whatever OS-version scheduling
        # behavior exists at the time, not assumptions made now.
        num_threads = max(1, joblib.cpu_count() - 1)

        for file_name in file_list:
            output_file = os.path.join(directory, blast_results_folder, file_name.replace(".junctions.fa", '.blast.txt'))
            fa_path = os.path.join(directory, blast_results_folder, file_name)
            total_sequences = self._count_fasta_sequences(fa_path)
            print(">>> Running BLAST search for file: %s (%d sequences)" % (file_name, total_sequences))
            blast_command_list = [os.path.join(blast_path, 'blastn' + suffix),
                                  '-query', os.path.join(directory, 'blast_results', file_name), '-db', db_path,
                                  '-task',  'blastn', '-dust', 'no', '-num_threads', str(num_threads),
                                  '-outfmt', '7', '-out', output_file, '-evalue', '0.2', '-max_target_seqs', '10']

            sys.stdout.flush()
            # A stale output_file from a previous failed attempt would
            # otherwise inflate the progress monitor's initial query count
            # until blastn's own fresh write overwrites it.
            if os.path.exists(output_file):
                os.remove(output_file)

            stop_event = threading.Event()
            progress_thread = threading.Thread(target=self._monitor_blast_progress,
                                               args=(output_file, total_sequences, stop_event),
                                               daemon=True)
            progress_thread.start()

            self.blast_pipe = subprocess.Popen(blast_command_list, shell=False, stderr=subprocess.PIPE)
            _, stderr_output = self.blast_pipe.communicate()
            stop_event.set()
            progress_thread.join(timeout=3.0)
            return_code = self.blast_pipe.returncode
            if return_code != 0 or not os.path.exists(output_file):
                # A failed blastn used to go unnoticed: nothing checked its
                # exit code (or captured its stderr), and
                # cleanup_junction_fa_files() deletes the source .fa
                # afterward regardless of whether a matching .blast.txt
                # actually got written - so a silent blastn crash meant
                # losing the .fa with no .blast.txt to show for it and no
                # error anywhere. Raising here stops the run before that
                # cleanup step and surfaces the failure, with blastn's own
                # error text, via junction_make_gui.py's excepthook
                # (junction_make_error.log in the work folder).
                raise RuntimeError(
                    "blastn failed (exit code %s) for %s. %s was not written.\nblastn stderr:\n%s" %
                    (return_code, file_name, output_file, stderr_output.decode('utf-8', errors='replace')))

    def generate_tabulated_blast_results(self, directory, blast_results_folder, blast_results_query_folder, gene_list_file, prefix):
        blast_list = [f for f in self.fileio.get_file_list(directory, blast_results_folder, ".txt") if f.startswith(prefix)]

        # Load the gene list once for the whole batch, rather than re-reading
        # and re-parsing the (often 100+MB) gene list file for every .blast.txt
        # file - it's the same list for every file in this batch.
        gene_list = self._get_accession_number_list(gene_list_file)

        for blasttxt in blast_list:
            print(">>> Parsing BLAST results file %s ..." % blasttxt)
            # The corresponding .fa file (still present - cleanup only runs
            # after every direction's parsing is done, see
            # cleanup_junction_fa_files()) gives an exact total for a
            # progress ETA, the same way blast_search() gets one.
            fa_path = os.path.join(directory, blast_results_folder,
                                   blasttxt.replace(".blast.txt", ".junctions.fa"))
            total_sequences = self._count_fasta_sequences(fa_path) if os.path.exists(fa_path) else None
            blast_dict, gene_dict = self._blast_parser(directory, blast_results_folder,
                                                       blasttxt, gene_list, total_sequences)
            for gene in list(blast_dict.keys()):
                if gene not in ['total', 'pos_que']:
                    stats = {'in_orf'  : 0, 'in_frame': 0, 'downstream': 0,
                             'upstream': 0, 'not_in_frame': 0,
                             'intron'  : 0, 'backwards': 0, 'frame_orf': 0, 'total': 0
                             }
                    for nm in blast_dict[gene].keys():
                        for j in blast_dict[gene][nm]:
                            j.ppm = j.count * 1000000.0 / blast_dict['total']
                            stats[j.frame] += 1
                            stats[j.orf] += 1
                            if j.frame_orf:
                                stats["frame_orf"] += 1
                            stats['total'] += 1
                    blast_dict[gene]['stats'] = stats

            blast_query_p = open(os.path.join(directory, blast_results_query_folder,
                                              blasttxt.replace(".blast.txt", ".bqp")), "wb")
            pickle.dump(blast_dict, blast_query_p)
            blast_query_p.close()

    def cleanup_junction_fa_files(self, directory, blast_results_folder, prefix):
        # Split out from generate_tabulated_blast_results() so BLAST, parsing,
        # and cleanup can each run as their own complete pass across every
        # enabled direction (5p and 3p), instead of blasting+parsing+cleaning
        # up one direction fully before starting the next.
        self.fileio.remove_file(directory, blast_results_folder,
                                [f for f in self.fileio.get_file_list(directory, blast_results_folder, ".fa") if f.startswith(prefix)])

    def _junctions_in_read(self, read, jseqs):
        match_index = -1
        junction_index = -1
        for i, j in enumerate(jseqs):
            if j in read:
                match_index = read.index(j)
                junction_index = i
        return (junction_index, match_index)

    def _search_for_junctions(self, input_filehandle, jseqs, exclusion_sequence, output_filehandle, f, input_file_size):
        hits_count = 0
        def check_matching_criteria(l, indexes, jseqs, read):
            value = 0
            if indexes[0] != -1:
                junction = jseqs[indexes[0]]
                downstream_rf = read[len(junction) + indexes[1] + (indexes[0] % 3) * 4:]
                if len(downstream_rf) > 25:
                    if exclusion_sequence not in downstream_rf or exclusion_sequence == '':
                        value = 1
                        protein_sequence = self.process.translate_orf(downstream_rf)
                        output_filehandle.write(" ".join(l[:4]) + " " + l[9] + " " + downstream_rf
                                                + " " + protein_sequence + "\n")
            return value

        reads_count = 0
        total_reads_count = 0
        while True:
            line = input_filehandle.readline()
            if not line:
                break
            line_split = line.strip().split()
            if line_split[0][0] != "@":
                total_reads_count += 1
                if line_split[2] == "*":
                    reads_count += 1
                    if reads_count % 5000 == 0:
                        sys.stdout.write('\rRead %.3f%% of file...' % (input_filehandle.tell() * 100.0 / input_file_size))
                        sys.stdout.flush()
                    sequence_read = line_split[9]
                    rev_sequence_read = self.process.reverse_complement(sequence_read)
                    fwd_indexes = self._junctions_in_read(sequence_read, jseqs)
                    hit = check_matching_criteria(line_split, fwd_indexes, jseqs, sequence_read)
                    if hit == 0:
                        rev_indexes = self._junctions_in_read(rev_sequence_read, jseqs)
                        hit = check_matching_criteria(line_split, rev_indexes, jseqs, rev_sequence_read)
                    hits_count += hit

        print(' ')  # terminate the \r-updating "Read X% of file..." line above
        print(">>> Found %d junctions in %d unmapped reads (%d total reads scanned)" %
             (hits_count, reads_count, total_reads_count))
        sys.stdout.flush()

    def _make_search_junctions(self, junctions_array):
        jseqs = []
        for junc in junctions_array:
            jseqs.append(junc[30:50])
            jseqs.append(junc[26:46])
            jseqs.append(junc[22:42])
        return jseqs

    def _get_accession_number_list(self, gene_list_file):
        fh = open(os.path.join('lists', gene_list_file), 'r')
        gene_list = {}
        for line in fh.readlines():
            split = line.split()
            gene_list[split[0]] = {'gene_name' : split[1],
                                   'orf_start' : int(split[6]) + 1,
                                   'orf_stop'  : int(split[7]),
                                   'mRNA'      : split[9],
                                   'intron'    : split[8],
                                   'chromosome': split[2]
                                   }
        return gene_list

    def _blast_parser(self, directory, infolder, fileName, gene_list, total_sequences=None):
        blast_results_handle = open(os.path.join(directory, infolder, fileName), 'r')
        blast_results_count = 0
        previous_bitscore = 0
        # Junctions are collected in a dict keyed by (position, query_start)
        # during parsing so a repeat hit is an O(1) lookup instead of an O(n)
        # scan of everything found so far for that gene/transcript. Converted
        # back to a plain list (in the same insertion order) before returning,
        # so the returned/pickled structure is unchanged.
        junction_lookup = {}
        gene_dict = {}
        collect_results = 'n'
        checkpoint = 10000
        last_checkpoint_count = 0
        last_checkpoint_time = time.time()
        printed_progress = False
        for line in blast_results_handle.readlines():
            line.strip()
            split = line.split()
            if "BLASTN" in line:
                blast_results_count += 1
                previous_bitscore = 0
                if total_sequences and blast_results_count - last_checkpoint_count >= checkpoint:
                    now = time.time()
                    elapsed = now - last_checkpoint_time
                    done_this_round = blast_results_count - last_checkpoint_count
                    rate = done_this_round / elapsed if elapsed > 0 else 0
                    remaining = max(0, total_sequences - blast_results_count)
                    eta_minutes = (remaining / rate / 60.0) if rate > 0 else 0
                    pct = blast_results_count * 100.0 / total_sequences
                    # Same \r-in-place convention as blastn's own progress
                    # line - see the matching comment in _monitor_blast_progress.
                    sys.stdout.write('\rParsing progress: %d / %d sequences (%.0f%%) - last %d took %.0fs, ~%.1f min remaining' %
                         (blast_results_count, total_sequences, pct, done_this_round, elapsed, eta_minutes))
                    sys.stdout.flush()
                    last_checkpoint_count = blast_results_count
                    last_checkpoint_time = now
                    printed_progress = True
            elif "hits" in line and int(split[1]) < 100:
                collect_results = 'y'
            elif split[0] != '#' and collect_results == 'y' and float(split[2]) > 98 and \
                            float(split[11]) > 50.0 and float(split[11]) > previous_bitscore:
                previous_bitscore = float(split[11]) * 0.98
                nm_number = split[1]
                gene_name = gene_list[nm_number]['gene_name']
                if gene_name not in gene_dict.keys():
                    gene_dict[gene_name] = [nm_number]
                else:
                    gene_dict[gene_name].append(nm_number)

                j = sts.jcnt()
                j.position = int(split[8])
                j.query_start = int(split[6])
                fudge_factor = j.query_start - 1
                frame = j.position - gene_list[nm_number]['orf_start'] - fudge_factor
                if frame % 3 == 0 or frame == 0:
                    j.frame = "in_frame"
                else:
                    j.frame = "not_in_frame"

                if gene_list[nm_number]['intron'] == "INTRON":
                    j.frame = "intron"

                if int(split[9]) - j.position < 0:
                    j.frame = "backwards"

                if j.position < gene_list[nm_number]['orf_start']:
                    j.orf = "upstream"
                elif j.position > gene_list[nm_number]['orf_stop']:
                    j.orf = "downstream"
                else:
                    j.orf = "in_orf"

                if j.frame == 'in_frame' and j.orf == 'in_orf':
                        j.frame_orf = True

                if gene_name not in junction_lookup:
                    junction_lookup[gene_name] = {}
                if nm_number not in junction_lookup[gene_name]:
                    junction_lookup[gene_name][nm_number] = {}

                key = (j.position, j.query_start)
                existing = junction_lookup[gene_name][nm_number].get(key)
                if existing is not None:
                    existing.count += 1
                else:
                    junction_lookup[gene_name][nm_number][key] = j
            else:
                collect_results = 'n'
        blast_results_handle.close()
        if printed_progress:
            print(' ')  # terminate the \r-updating progress line above

        results_dictionary = {}
        for gene_name, transcripts in junction_lookup.items():
            results_dictionary[gene_name] = {}
            for nm_number, hits_by_key in transcripts.items():
                results_dictionary[gene_name][nm_number] = list(hits_by_key.values())
        results_dictionary['total'] = blast_results_count
        return results_dictionary, gene_dict
