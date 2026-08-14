import sys
import os
import subprocess

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
     read for the junction/bait sequence(s) (pre-filled when you pick
     the library in step 1 - editable if needed).

     5' AND 3' JUNCTION SEQUENCE FIELDS: a library can have a 5'
     junction sequence, a 3' junction sequence, or both. Most libraries
     only have a 5' sequence on file (the 3' field is left blank); the
     two current TAB libraries (hORFeome_TAB/hg38 and SacCer3_TAB) have
     both. You can also paste in a 3' sequence by hand for a library
     that doesn't have one on file.

     BLAST 5'/3' JUNCTIONS CHECKBOXES: these control which of the two
     searches Junction Make actually runs, and are checked by default
     whenever their matching field has a sequence (a library with no 3'
     sequence on file has "Blast 3' Junctions" auto-unchecked, since
     there's nothing to search for - check the box yourself if you
     paste one in). You can run just the 5' search now and come back
     to run the 3' search on the same folder later (or vice versa) -
     each search only touches the files it created; re-running one
     direction never overwrites or deletes the other's output.

     Both searches use the exact same method - three overlapping
     20nt windows, BLAST filtering, .bqp generation, all described
     below - just against a different tag sequence. They run as two
     independent passes over the unmapped .sam files (once for 5', once
     for 3', if both are checked), and every file each pass produces is
     named with a "5p_" or "3p_" prefix so both directions' output can
     live side by side in the same folders without colliding. Blast
     Query, Read Depth, and Stat Maker don't need to know which
     directions you ran - they just work with whatever "5p_"/"3p_"
     files happen to exist.

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
     blast_results_query/5p_*.bqp and/or 3p_*.bqp - one file per
     .sam file per direction searched, a pickled per-gene,
     per-transcript map of every junction found (its position, frame,
     and whether it falls upstream/inside/downstream of the ORF).

  5. Blast Query lets you search a .bqp file by gene name (or NM_
     accession) and see/plot exactly where the junctions for that
     gene landed - position, frame, ORF status.

  6. Read Depth shows coverage depth across a gene's mRNA - search by
     gene name or accession, pick a dataset, and it counts, in
     windows across the transcript, how many reads (from the
     gene_count_indices data Gene Count already produced) overlap
     each window. A gene/dataset only works here once Gene Count has
     already processed that .sam file.

  7. Stat Maker (Analyze Data section - see its own Method Info button
     for the DESeq2 method and citation) compares a bait against a
     vector control using the gene_count_summary and .bqp files from
     steps 3-4, and produces one combined results table. Unlike the
     other Analyze Data tools, it needs no library or work folder
     selected here - it does its own folder/dataset selection - so
     its button is available as soon as DEEPN opens. It's also still
     available as its own separate app for anyone who just wants
     Stat Maker without the rest of DEEPN.

  8. FragFinder (Analyze Data section) focuses on one gene within one
     dataset at a time, lining up Read Depth, 5' junctions, and 3'
     junctions together instead of checking each tool separately.
     Enabled once the work folder has at least one .bqp file.

     GETTING STARTED: search for a gene (same gene name/accession
     search as Blast Query and Read Depth) and pick a transcript.
     Then, for any of the three panels (Read Depth, 5' Junctions, 3'
     Junctions), click that panel's "Load..." button and pick a file
     from the standard file dialog - a .sam file for Read Depth, a
     .bqp file for the junction panels. Picking a 5p_*.bqp/3p_*.bqp
     file (the naming Junction Make normally produces) for one panel
     automatically tries to fill in the other two by matching the
     base filename (swapping in the right prefix/suffix); if it
     can't find a match, it opens a file picker for you to locate it,
     and clicking Cancel there just leaves that panel empty. Older
     data with no 5p_/3p_ prefix at all works too - just pick each
     panel's file yourself, since there's no prefix for FragFinder to
     match a companion by. Only a 5' junction file is actually
     required; Read Depth and 3' junctions can be
     skipped for datasets that don't have them.

     Each panel's header shows which file is currently loaded (": "
     followed by the filename). The Read Depth panel has its own
     "Interval" spinbox (same as the standalone Read Depth tool) to
     adjust the window size used for the coverage calculation.

     All three plots share one position axis. The table on the
     right lists every 5' and 3' junction found for the current
     gene - click a row to mark that position with a dashed line on
     all three plots and enable "Extrapolate Junction".

     EXTRAPOLATE JUNCTION: the .bqp files only store an aggregated
     count per junction position, not the actual read sequences.
     For the junction you've selected, this re-scans that dataset's
     raw blast_results/*.blast.txt and junction_files/*.junctions.txt
     files (reapplying the same hit-quality filter Junction Make
     used) to find every individual read supporting that exact
     junction, and shows the longest one. Because it's scanning raw
     files that can be multiple gigabytes, this can take a couple of
     minutes on large datasets - the progress dialog updates
     continuously with how many MB have been scanned, so it stays
     visibly live even when it's slow.

     The result is two FASTA-formatted entries in the box on the
     right (Copy button available): first the junction/downstream
     sequence itself (gene name, 5' or 3', Forward or Backward),
     then the complete original sequencing read it came from (gene
     name, read identifier) - the underlying real sequencing read
     ID and sequence, not just the trimmed fragment that was BLASTed.

     The window can be resized: the junction table widens with the
     window, and the "Sequence of interest" box at the bottom
     absorbs any extra vertical space - the three plots stay a
     fixed size either way.

WHAT'S NEW IN DEEPN_26

Compared to the original 2018 build (Python 2 / PyQt4, Intel-only):

  - Ground-up port to Python 3 / PyQt5, running natively on Apple
    Silicon.

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

  - 3' junction support. Added a "3' Junction Sequence" field
    alongside the existing field (now labeled "5' Junction
    Sequence"), and "Blast 5' Junctions" / "Blast 3' Junctions"
    checkboxes to control which of the two Junction Make actually
    runs. See the Junction Make section above for the full behavior -
    both directions use the identical search/BLAST/.bqp method, run
    as independent passes, and can be run together or separately at
    different times without disturbing each other's output.

  - Junction Make status messages now say which direction (5' or 3')
    and which file is currently being searched, and print the actual
    tag windows being matched against, instead of a generic message.

  - New module: FragFinder. See step 8 above for full details - a
    single-dataset, single-gene view linking Read Depth, 5'
    junctions, and 3' junctions together, with a new "Extrapolate
    Junction" feature that finds the actual sequencing read behind a
    specific junction by re-scanning Junction Make's raw output.

  - Junction Make now finishes searching for junctions in every
    direction you've checked before BLASTing any of them, instead of
    fully finishing 5' (search through BLAST through cleanup) before
    3' searching even starts.

  - BLAST now leaves one CPU core free instead of claiming every one,
    so the rest of the machine (including DEEPN's own window) stays
    usable while it runs.

  - BLAST progress. blastn has no built-in progress option, but the
    growing results file itself shows how many sequences have been
    completed - so Junction Make now prints an update, with an
    estimated time remaining, every 10,000 completed sequences.

  - Aborting a running Gene Count/Junction Make/GC+JM (or quitting
    DEEPN entirely while one is running) now actually stops it,
    including any BLAST search in progress - previously the
    underlying processes, including blastn, kept running in the
    background with no way to stop them short of quitting them
    manually.

  - The live output box no longer greys out and becomes hard to read
    while Gene Count/Junction Make/GC+JM is running.

  - Blast Query, Read Depth, FragFinder, and Stat Maker can now be
    open at the same time as each other, instead of only one being
    allowed at once. They're still blocked while Gene Count/Junction
    Make/GC+JM is actually running, so BLAST always has the full
    machine to itself.

  - New: Stat Maker button (Analyze Data section). Needs no library
    or work folder selection - same as before, it does its own
    dataset selection - so it's available as soon as DEEPN opens.
    Still also available as its own separate app for anyone who just
    wants Stat Maker on its own.
"""

class vQlistWidgetItem(QtWidgets.QListWidgetItem):
    def __init__(self, value, data):
        super(vQlistWidgetItem, self).__init__('%s' % value)
        self.data = data

class DEEPN_Launcher(QtWidgets.QMainWindow, form_class):
    # Emitted from a background thread once check_path()'s os.path.exists()
    # checks finish - PyQt signals are safe to emit cross-thread, and
    # auto-marshal the connected slot call onto the receiver's own thread
    # (the main thread here), unlike calling widget methods directly from a
    # thread the way monitor_directory_for_changes() does.
    _existing_folders_found = QtCore.pyqtSignal(list)

    def __init__(self, *args):
        super(DEEPN_Launcher, self).__init__(*args)
        self.setupUi(self)
        self._existing_folders_found.connect(self._on_existing_folders_found)
        self._pending_check_done_callback = None
        self._caffeinate = None
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
        # status_layout (the live output box) is deliberately excluded here -
        # this list is only used to collect widgets that get disabled during
        # a run, and a disabled QPlainTextEdit renders with dimmed text,
        # making the running output hard to read right when it matters most.
        self._layouts = [self.analyze_data_layout, self.process_data_layout_1,
                         self.process_data_layout_2, self.process_data_layout_3,
                         self.process_data_layout_3p, self.process_data_layout_4]
        self.window = self.window()
        self.window.setGeometry(10, 30, self.width(), self.height())
        self.buttons = []
        self.thread = None
        # Message_Dialog sets Qt.WindowModal, which needs a parent to be
        # modal relative to - without one, macOS can't attach it as a
        # proper sheet/child window, which plausibly explains a hang seen
        # live: the dialog painted and was visible, but the process's own
        # event loop sat 100% idle in mach_msg2_trap (confirmed via
        # `sample`) - it was never actually receiving click events for
        # this window at all, only native window-manager-level actions
        # (the close button) worked.
        self.message = m.Message_Dialog(self)
        self.comment = m.Message_Dialog(self)
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

        # Viewer tools (Blast Query, Read Depth, FragFinder) each run in their
        # own QProcess and can be open concurrently with each other. They're
        # only blocked while a processing tool (Gene Count/Junction Make/
        # GC+JM, which share self.process below) is actually running, so
        # blastn always has the full machine to itself.
        self.viewer_buttons = [self.query_blast_btn, self.read_depth_btn, self.fragfinder_btn, self.stat_maker_btn]
        self._viewer_button_labels = {btn: btn.text() for btn in self.viewer_buttons}
        self._viewer_processes = {}

        # QProcess object for external app
        app.aboutToQuit.connect(self.on_about_to_quit)
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
        # Stat Maker needs no library/folder selection - it's usable
        # immediately, unlike every other Analyze Data button.
        self.stat_maker_btn.setEnabled(True)

        self.message.quit_btn.clicked.connect(self.message_quit_signal)
        self.message.continue_btn.clicked.connect(self.message_continue_signal)
        self.comment.continue_btn.clicked.connect(self.comment_continue_signal)
        self.comment.quit_btn.setEnabled(False)

    @QtCore.pyqtSlot()
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

    def _existing_folders(self, root, folders):
        # The actual os.path.exists() calls, pulled out of check_path() so
        # they can run on a background thread - on a slow/idle external or
        # network volume, each one is a real filesystem round-trip and a
        # handful of them run back-to-back was blocking the main thread
        # (and therefore the whole GUI, including the warning dialog this
        # feeds into) for several seconds, confirmed via a live test.
        return [d for d in folders if os.path.exists(os.path.join(root, d))]

    def check_path(self, existing_directories):
        # Must run on the main thread - self.message.exec_() starts its own
        # nested event loop. Takes an already-computed list; see
        # _check_path_async()/_on_existing_folders_found() for how callers
        # get here without blocking on the existence checks themselves.
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
            # raise_()/activateWindow() must run BEFORE exec_() - exec_()
            # blocks until the dialog closes, so calling them after (as this
            # used to) is a no-op. Showing the dialog synchronously from a
            # direct button-click handler let Qt/Cocoa activate the window
            # automatically as a natural continuation of that click; showing
            # it from _on_existing_folders_found() (a queued cross-thread
            # signal, not a direct user gesture) loses that automatic
            # activation - confirmed live: the dialog appeared with its
            # traffic-light buttons hollow/gray (inactive window), and only
            # the native close button (window-manager level) worked, not
            # Continue (needs the window to actually be key to receive it).
            self.message.raise_()
            self.message.activateWindow()
            self.message.exec_()

    def _check_path_async(self, folders, on_done):
        """Runs _existing_folders() on a background thread, then shows
        check_path()'s warning dialog (if needed) and calls on_done() -
        both back on the main thread, via _existing_folders_found - once
        the check completes."""
        self._pending_check_done_callback = on_done
        directory = self.directory
        def worker():
            self._existing_folders_found.emit(self._existing_folders(directory, folders))
        threading.Thread(target=worker, daemon=True).start()

    def _on_existing_folders_found(self, existing_directories):
        self.check_path(existing_directories)
        callback = self._pending_check_done_callback
        self._pending_check_done_callback = None
        if callback:
            callback()

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
            self.junction_sequence_3p_txt.setEnabled(True)
            self.exclude_sequence_txt.setEnabled(True)
            self.gene_count_btn.setEnabled(True)
            junction_ok = self.junction_prereqs_met()
            self.junction_make_btn.setEnabled(junction_ok)
            self.gene_count_junction_make_btn.setEnabled(junction_ok)
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

            self._update_viewer_buttons_enabled()

            if self.clicked_button != None:
                for btn in self.buttons1:
                    btn.setEnabled(False)
            else:
                for btn in self.buttons1:
                    btn.setEnabled(True)

            self.gene_count_btn.setText("Gene Count")
            self.junction_make_btn.setText("Junction Make")
            self.gene_count_junction_make_btn.setText("Gene Count + Junction Make")
            for btn in self.viewer_buttons:
                if btn not in self._viewer_processes:
                    btn.setText(self._viewer_button_labels[btn])
            time.sleep(0.1)


    def selection_changed(self):
        self.db_selection = self.db_list_wgt.item(self.db_list_wgt.currentRow())
        for row in self.g_db.select('*', id=str(self.db_selection.data)):
            self.junction_sequence_txt.setText(row[7])
            junction_3p = row[9] if len(row) > 9 and row[9] else ''
            self.junction_sequence_3p_txt.setText(junction_3p)
            self.blast_3p_box.setChecked(bool(junction_3p))
            self.gene_dictionary = str(row[4])
            self.chromosome_list = str(row[8])
            self.blast_db_name = str(row[6])
            self.gene_list_file = str(row[5])
        self.enable_select_folder()

    def enable_select_folder(self):
        for btn in self.buttons1:
            btn.setEnabled(True)

    def junction_prereqs_met(self):
        blast_5p = self.blast_5p_box.isChecked()
        blast_3p = self.blast_3p_box.isChecked()
        if not blast_5p and not blast_3p:
            return False
        if blast_5p and not str(self.junction_sequence_txt.text()).strip():
            return False
        if blast_3p and not str(self.junction_sequence_3p_txt.text()).strip():
            return False
        return True

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

        # process_started() (which calls this) sets self.proceed = 0, which
        # stops monitor_directory_for_changes() - the only other place this
        # gets called. Without this, the viewer buttons (Blast Query, Read
        # Depth, FragFinder, Stat Maker) stay frozen in whatever state they
        # were in the instant before processing started, i.e. enabled, for
        # the entire Gene Count/Junction Make/GC+JM run - even though that
        # run wants every core for blastn.
        self._update_viewer_buttons_enabled()

    def process_started(self):
        self.status_text.clear()
        self.proceed = 0
        self.disable_unused_buttons()
        # A long BLAST run can take hours; if the Mac sleeps overnight, an
        # open write handle to a file on an external/USB drive can be
        # silently killed by the kernel with no crash report and no Python
        # exception to catch - confirmed via `pmset -g log` after a real
        # overnight run died mid-write: repeated "Maintenance Sleep" cycles,
        # each one delayed ~1.1-1.2s waiting on IOUSBMassStorageInterfaceNub/
        # IOSCSIPeripheralDeviceType00 (the external drive's own USB
        # mass-storage driver) to acknowledge the transition.
        # `caffeinate -w <pid>` prevents idle/system sleep for exactly as
        # long as that specific process is alive, and exits on its own once
        # it is - doesn't touch the user's actual power settings.
        try:
            self._caffeinate = subprocess.Popen(
                ['caffeinate', '-i', '-s', '-w', str(self.process.processId())])
        except Exception:
            self._caffeinate = None

    def _stop_caffeinate(self):
        if self._caffeinate is not None:
            if self._caffeinate.poll() is None:
                self._caffeinate.terminate()
            self._caffeinate = None

    def process_finished(self):
        self._stop_caffeinate()
        self.quit = False
        self.status_bar.showMessage("%s process ended!" % self.clicked_button_text, 5000)
        self.clicked_button = None
        self.clicked_button_text = ''
        self.pid = None
        self.proceed = 1
        threading.Thread(target=self.monitor_directory_for_changes, daemon=True).start()

    def on_about_to_quit(self):
        # Without this, quitting DEEPN mid-run leaves whatever it launched
        # (Gene Count/Junction Make/GC+JM, or any open viewer tool) running
        # as an orphan - and for Junction Make specifically, its own blastn
        # child would be orphaned right along with it. terminate() sends
        # SIGTERM, which junction_make_gui.py now catches to clean up its
        # own blast subprocess before exiting.
        self._stop_caffeinate()
        if self.process.state() != QtCore.QProcess.NotRunning:
            self.process.terminate()
            self.process.waitForFinished(3000)
        for proc in list(self._viewer_processes.values()):
            if proc.state() != QtCore.QProcess.NotRunning:
                proc.terminate()
                proc.waitForFinished(3000)

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
            # Stop monitor_directory_for_changes() (a plain background thread
            # that mutates Qt widgets directly, not a QThread with proper
            # signal marshaling) before check_path()'s modal dialog runs its
            # own nested event loop below - otherwise the two race, which is
            # what was causing the warning dialog's Continue button to hang
            # unresponsive with a spinning cursor.
            self.proceed = 0
            self.disable_unused_buttons()
            def start():
                if self.quit == False:
                    self.status_bar.showMessage("Running %s ..." % self.clicked_button_text)
                    arguments = [self.directory, self.gene_dictionary, self.chromosome_list, str(self.combined)]
                    self.process.start(script_path('gene_count_gui.py'), arguments)
                else:
                    self.process_finished()
            if self.prompt == 2:
                self._check_path_async(['gene_count_summary', 'chromosome_files'], start)
            else:
                start()
        elif self.clicked_button == self.sender():
            self.process.terminate()
            self.kill_processes('gene_count_gui.py')

    @QtCore.pyqtSlot()
    def on_junction_make_btn_clicked(self):
        if self.clicked_button is None:
            self.clicked_button = self.sender()
            self.clicked_button_text = self.clicked_button.text()
            # See the matching comment in on_gene_count_btn_clicked.
            self.proceed = 0
            self.disable_unused_buttons()
            def start():
                if self.quit == False:
                    self.status_bar.showMessage("Running %s ..." % self.clicked_button_text)
                    arguments = [self.directory, str(self.junction_sequence_txt.text()),
                                 str(self.exclude_sequence_txt.text()), self.blast_db_name, self.gene_list_file, str(self.combined),
                                 str(self.junction_sequence_3p_txt.text()),
                                 str(int(self.blast_5p_box.isChecked())),
                                 str(int(self.blast_3p_box.isChecked()))]
                    self.process.start(script_path('junction_make_gui.py'), arguments)
                else:
                    self.process_finished()
            if self.prompt == 2:
                self._check_path_async(['junction_files', 'blast_results', 'blast_results_query'], start)
            else:
                start()
        elif self.clicked_button == self.sender():
            self.process.terminate()

    @QtCore.pyqtSlot()
    def on_gene_count_junction_make_btn_clicked(self):
        if self.clicked_button is None:
            self.clicked_button = self.sender()
            self.clicked_button_text = self.clicked_button.text()
            # See the matching comment in on_gene_count_btn_clicked.
            self.proceed = 0
            self.disable_unused_buttons()
            def start():
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
                                str(self.combined),
                                str(self.junction_sequence_3p_txt.text()),
                                str(int(self.blast_5p_box.isChecked())),
                                str(int(self.blast_3p_box.isChecked()))]
                    self.process.start(script_path('gc_jm.py'), arguments)
                else:
                    self.process_finished()
            if self.prompt == 2:
                self._check_path_async(['junction_files', 'blast_results',
                                        'blast_results_query', 'gene_count_summary',
                                        'chromosome_files'], start)
            else:
                start()
        elif self.clicked_button == self.sender():
            self.process.terminate()
            self.kill_processes('gene_count_gui.py')

    def _launch_viewer_tool(self, button, script, arguments):
        """Blast Query, Read Depth, and FragFinder each get their own QProcess
        (unlike Gene Count/Junction Make/GC+JM, which share self.process and
        stay mutually exclusive so blastn always has full use of the machine's
        cores). That means any number of viewer tools can be open together,
        as long as none of the processing tools are currently running - see
        _update_viewer_buttons_enabled()."""
        proc = self._viewer_processes.get(button)
        if proc is not None:
            proc.terminate()
            return
        proc = QtCore.QProcess(self)
        proc.finished.connect(lambda *a, b=button: self._on_viewer_process_finished(b))
        self._viewer_processes[button] = proc
        button.setText("Abort")
        self.status_bar.showMessage("Running %s ..." % self._viewer_button_labels[button])
        proc.start(script_path(script), arguments)

    def _on_viewer_process_finished(self, button):
        button.setText(self._viewer_button_labels[button])
        self.status_bar.showMessage("%s process ended!" % self._viewer_button_labels[button], 5000)
        self._viewer_processes.pop(button, None)

    def _viewer_prereq_met(self, button):
        if button is self.stat_maker_btn:
            # Stat Maker needs nothing from DEEPN's selected library/folder,
            # so it's the only one that doesn't depend on self.directory -
            # check it first, before self.directory is used below (None
            # until a folder is selected, e.g. at __init__ time).
            return True
        if self.directory is None:
            return False
        if button is self.query_blast_btn or button is self.fragfinder_btn:
            return os.path.exists(os.path.join(self.directory, 'blast_results_query')) and \
                len(self.fileio.get_file_list(self.directory, 'blast_results_query', '.bqp')) > 0
        if button is self.read_depth_btn:
            return os.path.exists(os.path.join(self.directory, 'gene_count_summary')) and \
                (len(self.fileio.get_file_list(self.directory, 'mapped_sam_files', '.sam')) > 0 or
                 len(self.fileio.get_file_list(self.directory, 'sam_files', '.sam')) > 0)
        return False

    def _update_viewer_buttons_enabled(self):
        processing_busy = self.clicked_button is not None
        for btn in self.viewer_buttons:
            if btn in self._viewer_processes:
                btn.setEnabled(True)
            else:
                btn.setEnabled((not processing_busy) and self._viewer_prereq_met(btn))

    @QtCore.pyqtSlot()
    def on_query_blast_btn_clicked(self):
        self._launch_viewer_tool(self.query_blast_btn, 'query_blast_gui.py',
                                 [self.directory, self.gene_list_file])

    @QtCore.pyqtSlot()
    def on_read_depth_btn_clicked(self):
        self._launch_viewer_tool(self.read_depth_btn, 'read_depth_gui_v2.py',
                                 [self.directory, self.gene_list_file, str(self.combined)])

    @QtCore.pyqtSlot()
    def on_fragfinder_btn_clicked(self):
        self._launch_viewer_tool(self.fragfinder_btn, 'fragfinder_gui.py',
                                 [self.directory, self.gene_list_file, str(self.combined)])

    @QtCore.pyqtSlot()
    def on_stat_maker_btn_clicked(self):
        self._launch_viewer_tool(self.stat_maker_btn, 'stat_maker_gui_v2.py', [])

form = DEEPN_Launcher()
form.show()
app.exec_()
app.deleteLater()
sys.exit()
