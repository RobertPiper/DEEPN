#!/usr/bin/env python3
import os
import sys

if '.app/Contents/MacOS' in sys.executable:
    os.chdir(os.path.join(os.path.dirname(sys.executable), '..', 'Resources'))

import signal
import time
import io
import traceback
import functions.fileio_gui as f
import functions.junctionf_gui as j
import functions.printio_gui as p


input_data_folder = 'unmapped_sam_files'     #Manage name of input folder here
junction_folder = 'junction_files'            #Manage name of junction reads output folder here
blast_results_folder = 'blast_results'            #Manage name of blast results output folder here
blast_results_query = 'blast_results_query'     #Manage name of blast results dictionary output folder here

main_directory = sys.argv[1]
junction_sequence = sys.argv[2].replace(" ", "").split(",")
exclusion_sequence = sys.argv[3].replace(" ", "")
blast_db_name = sys.argv[4]
gene_list_file = sys.argv[5]
combined = int(sys.argv[6])

junction_sequence_3p_raw = sys.argv[7].replace(" ", "") if len(sys.argv) > 7 else ''
junction_sequence_3p = junction_sequence_3p_raw.split(",") if junction_sequence_3p_raw else []
blast_5p = int(sys.argv[8]) if len(sys.argv) > 8 else 1
blast_3p = int(sys.argv[9]) if len(sys.argv) > 9 else 0

if combined == 1:
    input_data_folder = 'sam_files'

def initialize_folders(directory):
    fileio.create_new_folder(directory, junction_folder)
    # Creates the folder that will hold the Genecounts summaries
    fileio.create_new_folder(directory, blast_results_folder)
    # Creates the folder that will hold more granualar data on exon counts per chromosome
    fileio.create_new_folder(directory, blast_results_query)
    # printio.show_message("Files will be processed and saved to the following directories\n\n"
    #                      "1. Input Files : %s\n\n"
    #                      "2. Output Directory: %s\n\n" % (os.path.join(main_directory, input_folder),
    #                                                       os.path.join(main_directory, summary_folder)))


def excepthook(excType, excValue, tracebackobj):
    """
    Global function to catch unhandled exceptions.

    @param excType exception type
    @param excValue exception value
    @param tracebackobj traceback object
    """
    separator = '-' * 80
    logFile = main_directory + "/junction_make_error.log"
    notice = \
        """An unhandled exception occurred. Please report the problem\n""" \
        """using the error reporting dialog or via email to <%s>.\n""" \
        """A log has been written to "%s".\n\nError information:\n""" % \
        ("yourmail at server.com", "")

    timeString = time.strftime("%Y-%m-%d, %H:%M:%S")

    tbinfofile = io.StringIO()
    traceback.print_tb(tracebackobj, None, tbinfofile)
    tbinfofile.seek(0)
    tbinfo = tbinfofile.read()
    errmsg = '%s: \n%s' % (str(excType), str(excValue))
    sections = [separator, timeString, separator, errmsg, separator, tbinfo]
    msg = '\n'.join(sections)
    try:
        f = open(logFile, "w")
        f.write(msg)
        f.close()
    except IOError:
        pass

sys.excepthook = excepthook

if __name__ == '__main__':
    printio = p.printio()
    fileio = f.fileio()
    junctionf = j.junctionf(fileio, printio)
    sys.stdout.write("\n*** Junction Search ***\n\n")
    sys.stdout.flush()
    # signal.signal(signal.SIGTERM, junctionf.sigterm_handler)

    # printio.print_comment("Comment1")
    # fileio.input_data_check(main_directory, input_data_folder, '.sam',
    #                         [junction_folder, blast_results_folder, blast_results_query])
    initialize_folders(main_directory)

    if blast_5p:
        sys.stdout.write("\n>>> Running 5' junction search...\n")
        sys.stdout.flush()
        junctionf.junction_search(main_directory, junction_folder, input_data_folder, blast_results_folder,
                                 blast_results_query, junction_sequence, exclusion_sequence, '5p_')
        junctionf.blast_search(main_directory, blast_db_name, blast_results_folder, blast_results_query, '5p_')
        junctionf.generate_tabulated_blast_results(main_directory,
                                                   blast_results_folder,
                                                   blast_results_query,
                                                   gene_list_file, '5p_')

    if blast_3p and junction_sequence_3p:
        sys.stdout.write("\n>>> Running 3' junction search...\n")
        sys.stdout.flush()
        junctionf.junction_search(main_directory, junction_folder, input_data_folder, blast_results_folder,
                                 blast_results_query, junction_sequence_3p, exclusion_sequence, '3p_')
        junctionf.blast_search(main_directory, blast_db_name, blast_results_folder, blast_results_query, '3p_')
        junctionf.generate_tabulated_blast_results(main_directory,
                                                   blast_results_folder,
                                                   blast_results_query,
                                                   gene_list_file, '3p_')
