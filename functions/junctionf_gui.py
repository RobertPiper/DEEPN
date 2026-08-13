import os
import sys
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
        self._spin = None

    def sigterm_handler(self, _signo, _stack_frame):
        if self.blast_pipe:
            self.blast_pipe.terminate()

        if self._spin:
            self._spin.stop()

        if self.blast_pipe < 0:
            print(">>> Terminated Process (%d). Now Exiting Gracefully!" % self.blast_pipe)
            sys.stdout.flush()
        else:
            print(">>> All Process Terminated. Exiting Gracefully!")
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

        for file_name in file_list:
            output_file = os.path.join(directory, blast_results_folder, file_name.replace(".junctions.fa", '.blast.txt'))
            print(">>> Running BLAST search for file: " + file_name)
            blast_command_list = [os.path.join(blast_path, 'blastn' + suffix),
                                  '-query', os.path.join(directory, 'blast_results', file_name), '-db', db_path,
                                  '-task',  'blastn', '-dust', 'no', '-num_threads', str(joblib.cpu_count()),
                                  '-outfmt', '7', '-out', output_file, '-evalue', '0.2', '-max_target_seqs', '10']

            sys.stdout.flush()
            self.blast_pipe = subprocess.Popen(blast_command_list, shell=False, stderr=subprocess.PIPE)
            _, stderr_output = self.blast_pipe.communicate()
            return_code = self.blast_pipe.returncode
            if return_code != 0 or not os.path.exists(output_file):
                # A failed blastn used to go unnoticed: nothing checked its
                # exit code (or captured its stderr), and
                # generate_tabulated_blast_results() deletes the source .fa
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
            blast_dict, gene_dict = self._blast_parser(directory, blast_results_folder,
                                                       blasttxt, gene_list)
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
        while True:
            line = input_filehandle.readline()
            if not line:
                break
            line_split = line.strip().split()
            if line_split[0][0] != "@" and line_split[2] == "*":
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

    def _blast_parser(self, directory, infolder, fileName, gene_list):
        blast_results_handle = open(os.path.join(directory, infolder, fileName), 'r')
        blast_results_count = 0
        print_counter = 0
        previous_bitscore = 0
        # Junctions are collected in a dict keyed by (position, query_start)
        # during parsing so a repeat hit is an O(1) lookup instead of an O(n)
        # scan of everything found so far for that gene/transcript. Converted
        # back to a plain list (in the same insertion order) before returning,
        # so the returned/pickled structure is unchanged.
        junction_lookup = {}
        gene_dict = {}
        collect_results = 'n'
        for line in blast_results_handle.readlines():
            line.strip()
            split = line.split()
            if "BLASTN" in line:
                blast_results_count += 1
                print_counter += 1
                previous_bitscore = 0
                if print_counter == 90000:
                    sys.stdout.write('.')
                    print_counter = 0
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

        results_dictionary = {}
        for gene_name, transcripts in junction_lookup.items():
            results_dictionary[gene_name] = {}
            for nm_number, hits_by_key in transcripts.items():
                results_dictionary[gene_name][nm_number] = list(hits_by_key.values())
        results_dictionary['total'] = blast_results_count
        return results_dictionary, gene_dict
