import sys
import os

if '.app/Contents/MacOS' in sys.executable:
    os.chdir(os.path.join(os.path.dirname(sys.executable), '..', 'Resources'))

import re
import time
import threading
import itertools
import functions.db as db
import functions.fileio_gui as f
import functions.printio_gui as p
import functions.message as m

from PyQt5 import QtCore, QtGui, QtWidgets, uic

app = QtWidgets.QApplication(sys.argv)
ui_path = os.path.join('ui', 'DEEPN.ui')
if sys.platform == 'win32':
    ui_path = os.path.join('ui', 'Windows', 'DEEPN.ui')
form_class, base_class = uic.loadUiType(ui_path)

# In dev mode (running from source rather than a packaged .app bundle), the
# launcher spawns the other tools as plain Python scripts from this directory.
# In a packaged .app, they are built as separate executables (via py2app's
# extra_scripts) living alongside the main executable in Contents/MacOS/.
FROZEN = '.app/Contents/MacOS' in sys.executable
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MACOS_DIR = os.path.dirname(sys.executable)

def script_path(name):
    if FROZEN:
        return os.path.join(MACOS_DIR, os.path.splitext(name)[0])
    return os.path.join(SCRIPT_DIR, name)

INFO_TEXT = """DEEPN - LAUNCHER INFO

DEEPN is a workflow for analyzing yeast two-hybrid (Y2H) next-generation
sequencing data.

HOW TO USE THIS

  1. Select a Y2H library/database from the list on the left. This
     picks which genome you're working with (e.g. mouse mm10, human
     hg38, yeast sacCer3) and loads everything tied to it: the gene
     dictionary, chromosome list, BLAST reference database, gene list,
     and the junction/bait sequence used to find junction reads.

  2. Click "Select Folder" and choose (or create) your work folder.
     DEEPN creates three subfolders in it if they don't exist:
       mapped_sam_files    - .sam files of reads that mapped to the
                              genome (input to Gene Count)
       unmapped_sam_files  - .sam files of reads that did NOT map
                              (input to Junction Make)
       sam_files           - if you're running combined mode, both
                              mapped and unmapped reads together go
                              here instead of the two folders above
     Put your sequencing output .sam files into the appropriate
     folder(s) before running anything.

  3. Gene Count reads the mapped .sam files and counts how many reads
     land within each gene's annotated exons. It writes:
       gene_count_summary/*.csv   - PPM (parts-per-million) per gene,
                                    used by Read Depth's gene search,
                                    Stat Maker, and library QC
       chromosome_files/*.txt     - per-exon detail, including raw
                                    read counts
       gene_count_indices/*/*.bin - per-chromosome read position data,
                                    used by Read Depth to calculate
                                    coverage

  4. Junction Make reads the unmapped .sam files and searches each
     read for the junction/bait sequence (pre-filled when you pick
     the library in step 1 - editable if needed).

     JUNCTION SEARCH TOLERANCE: rather than searching for one fixed
     20nt window right at the junction boundary, three overlapping
     20nt windows are searched, each shifted 4nt earlier than the
     last. This gives tolerance for a sequencing error landing right
     at one particular window's edge - a real junction can still be
     found via one of the other two offsets, and the downstream
     reading frame is adjusted to match whichever offset actually hit.

     Where a match is found, the sequence downstream of the junction
     is pulled out as a candidate coding fragment and written to a
     FASTA file. THIS IS WHERE BLASTN RUNS: each of those fragments is
     searched with blastn against the genome's reference transcriptome
     database to identify which gene/transcript it actually came from.

     BLAST HIT FILTER: not every blastn match is kept. A hit is only
     counted if: the read had fewer than 100 total BLAST hits (reads
     matching too many places are considered too ambiguous to use),
     percent identity is above 98%, bitscore is above 50, and bitscore
     is within 2% of the best bitscore seen so far for that read (so
     near-tied top hits are kept, not just a single best match, but
     weak/spurious hits are dropped).

     The surviving hits are parsed and written to
     blast_results_query/*.bqp - one file per .sam file, a pickled
     per-gene, per-transcript map of every junction found (its
     position, frame, and whether it falls upstream/inside/downstream
     of the ORF).

  5. Blast Query lets you search a .bqp file by gene name (or NM_
     accession) and see/plot exactly where the junctions for that
     gene landed - position, frame, ORF status.

  6. Read Depth shows coverage depth across a gene's mRNA - search by
     gene name or accession, pick a dataset, and it counts, in
     windows across the transcript, how many reads (from the
     gene_count_indices data Gene Count already produced) overlap
     each window. A gene/dataset only works here once Gene Count has
     already processed that .sam file.

  7. Stat Maker (launched separately - see its own Method Info button)
     compares a bait against a vector control using the
     gene_count_summary and .bqp files from steps 3-4, and produces
     one combined results table.

WHAT'S NEW IN THIS VERSION (DEEPN_26v1)

This is a ground-up port from the original Python 2 / PyQt4 build (last
updated 2018, Intel-only) to Python 3 / PyQt5, running natively on Apple
Silicon. Beyond the language/framework port itself:

  - BLAST now runs natively. The original bundled an Intel-only BLAST
    2.2.30 binary, which needed Rosetta translation to run at all on
    Apple Silicon. This version uses a native arm64 BLAST+ 2.17.0
    installation, with search databases rebuilt from the reference gene
    lists - no translation layer involved anywhere in Junction Make.

  - Read Depth bug fix. The reverse-complement routine used to
    calculate coverage crashed on lowercase (soft-masked) bases in
    reference mRNA sequences - a real bug in the original code, just
    never triggered until tested against real hg38 data. Fixed to be
    case-insensitive.

  - Plotting bug fix. A newer version of the plotting library
    (pyqtgraph) requires hex color codes to start with '#'; the
    original bare hex strings crashed both Blast Query's and Read
    Depth's plots. Fixed.

  - Read Depth: gene-name search. You can now search by gene name and
    pick from its list of NM_ transcript variants, instead of only
    being able to search by exact NM_ accession number.

  - Stat Maker rebuilt. JAGS-based Bayesian statistics have been
    retired entirely. Statistics are now computed with DESeq2 (see
    Stat Maker's own Method Info for the full method and citation),
    and junction-data collation from BLAST results was rebuilt as an
    independent step.

  - Recovered reference data. The gene list files (lists/*.prn)
    needed by Blast Query, Read Depth, and Junction Make were missing
    from the source repository and have been restored from the
    original DEEPN 3.7 build.

WHAT'S NEW IN v3

  - Faster BLAST-results-to-.bqp conversion (Junction Make). Two
    inefficiencies in generating blast_results_query/*.bqp files from
    blastn's output were fixed, with no change to the output format
    or the filtering criteria (see the Junction Make section above
    for what those filters are - unchanged):

      1. The gene list (lists/*.prn, often 100+ MB) used to be
         re-read and re-parsed from disk for every single .blast.txt
         file being processed. It's now loaded once per batch and
         reused, since it's the same list for every file.

      2. Duplicate-junction detection used to scan every junction
         already found for a gene/transcript, for every new BLAST hit
         line - meaning the check got slower the more hits a gene had
         already accumulated. It's now a direct lookup instead, with
         no change to which junctions get counted or how they're
         counted (still by exact position + query start, still a
         simple +1 to an existing count on a repeat).
"""

class vQlistWidgetItem(QtWidgets.QListWidgetItem):
    def __init__(self, value, data):
        super(vQlistWidgetItem, self).__init__('%s' % value)
        self.data = data

class DEEPN_Launcher(QtWidgets.QMainWindow, form_class):
    def __init__(self, *args):
        super(DEEPN_Launcher, self).__init__(*args)
        self.setupUi(self)
        self.proceed = 1
        self.prompt = 2
        self.fileio = f.fileio()
        self.printio = p.printio()
        self.directory = None
        self.clicked_button = None
        self.clicked_button_text = None
        self.db_selection = None
        self.combined = 0
        self.blast_db_name = ''
        self.gene_dictionary = ''
        self.gene_list_file = ''
        self.chromosome_list = ''
        self.start_match = re.compile(r'^[>>>|\t|***]')
        self.bar = itertools.cycle(['/', '-', '\\'])
        self._layouts = [self.analyze_data_layout, self.process_data_layout_1,
                         self.process_data_layout_2, self.process_data_layout_3,
                         self.process_data_layout_4, self.status_layout]
        self.window = self.window()
        self.window.setGeometry(10, 30, self.width(), self.height())
        self.buttons = []
        self.thread = None
        self.message = m.Message_Dialog()
        self.comment = m.Message_Dialog()
        self.quit = False
        for layout in self._layouts:
            widgets = (layout.itemAt(i).widget() for i in range(layout.count()))
            for btn in widgets:
                if btn != None:
                    self.buttons.append(btn)

        self.buttons1 = []
        widgets = (self.layout1.itemAt(i).widget() for i in range(self.layout1.count()))
        for btn in widgets:
            if btn != None:
                self.buttons1.append(btn)

        # Checkbox
        self.prompt_box.stateChanged.connect(self.on_prompt_box_stateChanged)
        self.info_btn.clicked.connect(self.on_info_btn_clicked)

        # QProcess object for external app
        self.process = QtCore.QProcess(self)
        self.process.readyReadStandardOutput.connect(self.stdout_ready)
        # self.process.readyReadStandardError.connect(self.stderr_ready)
        self.process.started.connect(self.process_started)
        self.process.finished.connect(self.process_finished)

        # QThread for processes
        self.thread = QtCore.QThread(self)

        # Connect database
        self.g_db = db.genome_db(os.path.join('DEEPN_db.sqlite3'))
        for row in self.g_db.select_all('*'):
            self.db_list_wgt.addItem(vQlistWidgetItem(row[1], row[0]))

        # Binding Selection Changed Signal
        self.db_list_wgt.itemSelectionChanged.connect(self.selection_changed)

        # Disable All buttons
        self.disable_unused_buttons()
        for btn in self.buttons1:
            btn.setEnabled(False)

        self.message.quit_btn.clicked.connect(self.message_quit_signal)
        self.message.continue_btn.clicked.connect(self.message_continue_signal)
        self.comment.continue_btn.clicked.connect(self.comment_continue_signal)
        self.comment.quit_btn.setEnabled(False)

    def on_info_btn_clicked(self):
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("DEEPN Info")
        dialog.resize(800, 700)
        layout = QtWidgets.QVBoxLayout(dialog)
        text = QtWidgets.QPlainTextEdit(dialog)
        text.setReadOnly(True)
        text.setFont(QtGui.QFont("Menlo", 12))
        text.setPlainText(INFO_TEXT)
        layout.addWidget(text)
        close_btn = QtWidgets.QPushButton("Close", dialog)
        close_btn.clicked.connect(dialog.accept)
        layout.addWidget(close_btn)
        dialog.exec_()

    def message_quit_signal(self):
        self.quit = True
        self.message.hide()

    def message_continue_signal(self):
        self.quit = False
        self.message.hide()

    def comment_continue_signal(self):
        self.comment.hide()

    def check_and_create_initial_folders(self):
        if os.path.exists(os.path.join(self.directory, 'mapped_sam_files')):
            pass
        else:
            os.makedirs(os.path.join(self.directory, 'mapped_sam_files'))

        if os.path.exists(os.path.join(self.directory, 'unmapped_sam_files')):
            pass
        else:
            os.makedirs(os.path.join(self.directory, 'unmapped_sam_files'))

        if os.path.exists(os.path.join(self.directory, 'sam_files')):
            pass
        else:
            os.makedirs(os.path.join(self.directory, 'sam_files'))

    def check_path(self, root, folders):
        existing_directories = []
        for dir in folders:
            if os.path.exists(os.path.join(root, dir)):
                existing_directories.append(dir)
        if len(existing_directories) > 0:
            message = "One or more of the following folders is already present:<br>"
            count = 1
            for d in existing_directories:
                message += "<br>%d.<b>%s</b><br>" % (count, d)
                count += 1
            message += "<br><br>These folders and/or their contexts can be moved to avoid risk of being over-written"
            self.message.windowTitle('Warning!')
            self.message.showMessage('<b>WARNING</b><br><br>%s<br>' % message)
            self.message.continue_btn.setEnabled(True)
            self.message.exec_()
            self.message.activateWindow()

    def print_comment(self, tag):
        text = self.printio.get_text_block(tag)
        self.comment.windowTitle('Message')
        self.comment.showMessage('%s' % text)
        self.comment.continue_btn.setEnabled(True)
        self.comment.show()

    def get_main_directory(self):
        self.directory = str(QtWidgets.QFileDialog.getExistingDirectory(QtWidgets.QFileDialog(), "Locate Work Folder",
                                                                    os.path.expanduser("~"),
                                                                    QtWidgets.QFileDialog.ShowDirsOnly))
        if self.directory == None:
            self.get_main_directory()
        else:
            self.check_and_create_initial_folders()
            self.disable_unused_buttons()
            self.proceed = 1
            threading.Thread(target=self.monitor_directory_for_changes, daemon=True).start()
            pass

    def monitor_directory_for_changes(self):
        def enable_buttons():
            self.junction_sequence_txt.setEnabled(True)
            self.exclude_sequence_txt.setEnabled(True)
            self.gene_count_btn.setEnabled(True)
            self.junction_make_btn.setEnabled(True)
            self.gene_count_junction_make_btn.setEnabled(True)
            self.db_list_wgt.setEnabled(True)
            self.status_txt.setText(self.directory)

        while self.proceed == 1:
            if len(self.fileio.get_file_list(self.directory, 'sam_files', '.sam')) == 0:
                if len(self.fileio.get_file_list(self.directory, 'mapped_sam_files', '.sam')) <= 0 or \
                        len(self.fileio.get_file_list(self.directory, 'unmapped_sam_files', '.sam')) <= 0:
                    self.gene_count_btn.setEnabled(False)
                    self.junction_make_btn.setEnabled(False)
                    self.gene_count_junction_make_btn.setEnabled(False)
                    self.status_txt.setText('Waiting to load .sam files into selected directory...')
                else:
                    self.combined = 0
                    enable_buttons()
            else:
                self.combined = 1
                enable_buttons()

            if os.path.exists(os.path.join(self.directory, 'blast_results_query')):
                if len(self.fileio.get_file_list(self.directory, 'blast_results_query', '.bqp')) <= 0:
                    self.query_blast_btn.setEnabled(False)
                else:
                    self.query_blast_btn.setEnabled(True)
            else:
                self.query_blast_btn.setEnabled(False)


            if os.path.exists(os.path.join(self.directory, 'gene_count_summary')) and \
                    (len(self.fileio.get_file_list(self.directory, 'mapped_sam_files', '.sam')) > 0 or
                     len(self.fileio.get_file_list(self.directory, 'sam_files', '.sam')) > 0):
                self.read_depth_btn.setEnabled(True)
            else:
                self.read_depth_btn.setEnabled(False)

            if self.clicked_button != None:
                for btn in self.buttons1:
                    btn.setEnabled(False)
            else:
                for btn in self.buttons1:
                    btn.setEnabled(True)

            self.gene_count_btn.setText("Gene Count")
            self.junction_make_btn.setText("Junction Make")
            self.gene_count_junction_make_btn.setText("Gene Count + Junction Make")
            self.query_blast_btn.setText("Blast Query")
            self.read_depth_btn.setText("Read Depth")
            time.sleep(0.1)


    def selection_changed(self):
        self.db_selection = self.db_list_wgt.item(self.db_list_wgt.currentRow())
        for row in self.g_db.select('*', id=str(self.db_selection.data)):
            self.junction_sequence_txt.setText(row[7])
            self.gene_dictionary = str(row[4])
            self.chromosome_list = str(row[8])
            self.blast_db_name = str(row[6])
            self.gene_list_file = str(row[5])
        self.enable_select_folder()

    def enable_select_folder(self):
        for btn in self.buttons1:
            btn.setEnabled(True)

    def stdout_ready(self):
        text = bytes(self.process.readAllStandardOutput()).decode('utf-8', errors='replace').strip()
        self.append(text)

    def stderr_ready(self):
        text = bytes(self.process.readAllStandardError()).decode('utf-8', errors='replace')
        self.append(text)

    def append(self, text):
        self.status_bar.showMessage("Running %s script...  %s" % (self.clicked_button_text, next(self.bar)))
        cursor = self.status_text.textCursor()
        if self.start_match.match(text):
            cursor.select(QtGui.QTextCursor.LineUnderCursor)
            cursor.removeSelectedText()
            cursor.deleteChar()
            cursor.movePosition(QtGui.QTextCursor.EndOfLine)
            cursor.insertText(text + "\n\n")
        else:
            cursor.select(QtGui.QTextCursor.LineUnderCursor)
            cursor.removeSelectedText()
            cursor.deleteChar()
            cursor.insertText(text)
        self.status_text.ensureCursorVisible()

    def disable_unused_buttons(self):
        for btn in self.buttons:
            if btn is not self.clicked_button:
                btn.setEnabled(False)
            else:
                btn.setText("Abort")

        for btn in self.buttons1:
            btn.setEnabled(False)

        if self.clicked_button:
            self.db_list_wgt.setEnabled(False)

    def process_started(self):
        self.status_text.clear()
        self.proceed = 0
        self.disable_unused_buttons()

    def process_finished(self):
        self.quit = False
        self.status_bar.showMessage("%s process ended!" % self.clicked_button_text, 5000)
        self.clicked_button = None
        self.clicked_button_text = ''
        self.pid = None
        self.proceed = 1
        threading.Thread(target=self.monitor_directory_for_changes, daemon=True).start()

    def kill_processes(self, name):
        time.sleep(1)
        os.system('kill $(ps aux | awk \'/' + name + '/ {print $2}\')')

    @QtCore.pyqtSlot()
    def on_prompt_box_stateChanged(self):
        self.prompt = self.prompt_box.checkState()

    @QtCore.pyqtSlot()
    def on_select_folder_btn_clicked(self):
        self.get_main_directory()

    @QtCore.pyqtSlot()
    def on_gene_count_btn_clicked(self):
        print(self.prompt)
        if self.clicked_button is None:
            self.clicked_button = self.sender()
            self.clicked_button_text = self.clicked_button.text()
            if self.prompt == 2:
                self.check_path(self.directory, ['gene_count_summary', 'chromosome_files'])
            if self.quit == False:
                self.status_bar.showMessage("Running %s ..." % self.clicked_button_text)
                arguments = [self.directory, self.gene_dictionary, self.chromosome_list, str(self.combined)]
                self.process.start(script_path('gene_count_gui.py'), arguments)
            else:
                self.process_finished()
        elif self.clicked_button == self.sender():
            self.process.terminate()
            self.kill_processes('gene_count_gui.py')

    @QtCore.pyqtSlot()
    def on_junction_make_btn_clicked(self):
        if self.clicked_button is None:
            self.clicked_button = self.sender()
            self.clicked_button_text = self.clicked_button.text()
            if self.prompt == 2:
                self.check_path(self.directory, ['junction_files', 'blast_results', 'blast_results_query'])
            if self.quit == False:
                self.status_bar.showMessage("Running %s ..." % self.clicked_button_text)
                arguments = [self.directory, str(self.junction_sequence_txt.text()),
                             str(self.exclude_sequence_txt.text()), self.blast_db_name, self.gene_list_file, str(self.combined)]
                self.process.start(script_path('junction_make_gui.py'), arguments)
            else:
                self.process_finished()
        elif self.clicked_button == self.sender():
            self.process.terminate()

    @QtCore.pyqtSlot()
    def on_gene_count_junction_make_btn_clicked(self):
        if self.clicked_button is None:
            self.clicked_button = self.sender()
            self.clicked_button_text = self.clicked_button.text()
            if self.prompt == 2:
                self.check_path(self.directory, ['junction_files', 'blast_results',
                                                 'blast_results_query', 'gene_count_summary',
                                                 'chromosome_files'])
            if self.quit == False:
                self.status_bar.showMessage("Running %s ..." % self.clicked_button_text)
                arguments = [self.directory,
                            self.gene_dictionary, self.chromosome_list,
                            str(self.junction_sequence_txt.text()),
                            str(self.exclude_sequence_txt.text()),
                            self.blast_db_name,
                            script_path('gene_count_gui.py'),
                            script_path('junction_make_gui.py'),
                            self.gene_list_file,
                            str(self.combined)]

                self.process.start(script_path('gc_jm.py'), arguments)
            else:
                self.process_finished()
        elif self.clicked_button == self.sender():
            self.process.terminate()
            self.kill_processes('gene_count_gui.py')

    @QtCore.pyqtSlot()
    def on_query_blast_btn_clicked(self):
        if self.clicked_button is None:
            self.clicked_button = self.sender()
            self.clicked_button_text = self.clicked_button.text()
            self.status_bar.showMessage("Running %s ..." % self.clicked_button_text)
            arguments = [self.directory, self.gene_list_file]
            self.process.start(script_path('query_blast_gui.py'), arguments)
        elif self.clicked_button == self.sender():
            self.process.terminate()

    @QtCore.pyqtSlot()
    def on_read_depth_btn_clicked(self):
        if self.clicked_button is None:
            self.clicked_button = self.sender()
            self.clicked_button_text = self.clicked_button.text()
            self.status_bar.showMessage("Running %s ..." % self.clicked_button_text)
            self.process.start(script_path('read_depth_gui_v2.py'), [self.directory,
                                                self.gene_list_file, str(self.combined)])
        elif self.clicked_button == self.sender():
            self.process.terminate()

form = DEEPN_Launcher()
form.show()
app.exec_()
app.deleteLater()
sys.exit()
