#!/usr/bin/env python3
import sys
import signal
import subprocess

main_directory = sys.argv[1]
mapster_exe = sys.argv[2]
mapster_config_path = sys.argv[3]
gene_dictionary = sys.argv[4]
chromosomes_list_name = sys.argv[5]

junction_sequence = sys.argv[6]
exclusion_sequence = sys.argv[7]
blast_db_name = sys.argv[8]

gene_count_exe = sys.argv[9]
junction_make_exe = sys.argv[10]

gene_list_file = sys.argv[11]
combined = sys.argv[12]

junction_sequence_3p = sys.argv[13] if len(sys.argv) > 13 else ''
blast_5p = sys.argv[14] if len(sys.argv) > 14 else '1'
blast_3p = sys.argv[15] if len(sys.argv) > 15 else '0'

current_pipe = None

def sigterm_handler(_signo, _stack_frame):
    # Same reasoning as gc_jm.py's own handler: without this, SIGTERM (from
    # DEEPN's Abort button, or DEEPN quitting mid-run) kills this process
    # immediately with no chance to stop whichever of MAPster/Gene Count/
    # Junction Make is currently running as our own child - it would just
    # keep running, orphaned (and for MAPster specifically, its own
    # in-progress hisat2 alignment would be orphaned right along with it).
    # Junction Make itself already forwards SIGTERM to clean up its own
    # blastn child in turn.
    if current_pipe is not None and current_pipe.poll() is None:
        current_pipe.terminate()
    sys.exit()

signal.signal(signal.SIGTERM, sigterm_handler)

# MAPster runs its own GUI for picking input files and naming outputs (not
# something that can be predetermined - see the DEEPN-config loading in
# LOCAL_MAPster's mainwindow.cpp for what *is* pre-set: genome selection
# and output folder). It quits itself once its job queue empties, so this
# wait() unblocks exactly when a normal user-driven run there is done -
# same completion signal as every other step here, no new IPC needed.
current_pipe = subprocess.Popen([mapster_exe, mapster_config_path], shell=False)
current_pipe.wait()

current_pipe = subprocess.Popen([gene_count_exe, main_directory, gene_dictionary, chromosomes_list_name, combined], shell=False)
current_pipe.wait()

current_pipe = subprocess.Popen([junction_make_exe, main_directory, junction_sequence, exclusion_sequence, blast_db_name, gene_list_file, combined,
                                 junction_sequence_3p, blast_5p, blast_3p], shell=False)
current_pipe.wait()
