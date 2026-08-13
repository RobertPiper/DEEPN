#!/usr/bin/env python3
import sys
import signal
import subprocess

main_directory = sys.argv[1]
gene_dictionary = sys.argv[2]
chromosomes_list_name = sys.argv[3]

junction_sequence = sys.argv[4]
exclusion_sequence = sys.argv[5]
blast_db_name = sys.argv[6]

gene_count_exe = sys.argv[7]
junction_make_exe = sys.argv[8]

gene_list_file = sys.argv[9]
combined = sys.argv[10]

junction_sequence_3p = sys.argv[11] if len(sys.argv) > 11 else ''
blast_5p = sys.argv[12] if len(sys.argv) > 12 else '1'
blast_3p = sys.argv[13] if len(sys.argv) > 13 else '0'

current_pipe = None

def sigterm_handler(_signo, _stack_frame):
    # Without this, SIGTERM (from DEEPN's Abort button, or DEEPN quitting
    # mid-run) kills this process immediately with no chance to stop
    # whichever of Gene Count/Junction Make is currently running as our
    # own child - it would just keep running, orphaned. Junction Make
    # itself now forwards SIGTERM the same way to clean up its own blastn
    # child in turn.
    if current_pipe is not None and current_pipe.poll() is None:
        current_pipe.terminate()
    sys.exit()

signal.signal(signal.SIGTERM, sigterm_handler)

current_pipe = subprocess.Popen([gene_count_exe, main_directory, gene_dictionary, chromosomes_list_name, combined], shell=False)
current_pipe.wait()

current_pipe = subprocess.Popen([junction_make_exe, main_directory, junction_sequence, exclusion_sequence, blast_db_name, gene_list_file, combined,
                                 junction_sequence_3p, blast_5p, blast_3p], shell=False)
current_pipe.wait()
