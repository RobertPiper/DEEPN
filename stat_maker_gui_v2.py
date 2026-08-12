#!/usr/bin/env python3
import os
import sys

if '.app/Contents/MacOS' in sys.executable:
    os.chdir(os.path.join(os.path.dirname(sys.executable), '..', 'Resources'))

import re
import time
import glob
from datetime import datetime

from PyQt5 import QtCore, QtGui, QtWidgets, uic
import xlsxwriter as xls

import functions.fileio_gui as f
import functions.stat_collation as sc
import DragDropListView

app = QtWidgets.QApplication(sys.argv)
form_class, base_class = uic.loadUiType(os.path.join('ui', 'Stat_Maker_v2.ui'))

R_SCRIPT_PATH = os.path.join('r_scripts', 'run_Y2H_enrichement_stats.R')

METHOD_INFO_TEXT = """STAT MAKER v2 - METHOD OVERVIEW

WHAT THIS PROGRAM DOES

Stat Maker compares a Y2H bait (Selected vs Non-Selected) against a Vector
control (2 replicates, Selected vs Non-Selected each), then combines two
independent results into one table, one row per gene:

  1. STATISTICS (via DESeq2, see R code below)
     For each gene: pvalue, log2FoldChange, and a rank-based Enrichment_score,
     computed separately for Bait (Selected vs Non-Selected) and for Vector
     (Selected vs Non-Selected).

  2. JUNCTION COLLATION (from BLAST results, .bqp files)
     For each gene, in each of the 6 input datasets, the percentage of
     junctions that are:
       inframe_inorf - in-frame AND within the annotated ORF
       upstream      - start upstream of the ORF
       in_orf        - within the ORF (regardless of frame)
       downstream    - start downstream of the ORF
       in_frame      - in-frame (regardless of ORF position)
       backwards     - mapped in the reverse orientation
       intron        - ORF contains an intron (yeast genes only; rare)

INPUT

  Gene-level counts come from DEEPN's gene_count_summary CSVs, which report
  PPM (parts per million of total mapped reads), not raw read counts. Genes
  annotated on multiple scaffolds/alt-contigs are summed to one PPM total
  per gene before anything else happens.

STATISTICS METHOD (DESeq2)

  This DESeq2-based approach is adapted from the method described in:

    Velasquez-Zapata V, Elmore JM, Banerjee S, Dorman KS, Wise RP (2021)
    Next-generation yeast-two-hybrid analysis with Y2H-SCORES identifies
    novel interactors of the MLA immune receptor. PLOS Computational
    Biology 17(4): e1008890. https://doi.org/10.1371/journal.pcbi.1008890

  That method has been streamlined here to fit within the Stat Maker
  workflow (single collated input file, PPM-based counts, single bait vs.
  vector comparison).

  The count matrix handed to DESeq2 is PPM, not raw counts. Because PPM is
  already normalized for sequencing depth, DESeq2's own size-factor
  normalization is bypassed (normalized=TRUE in the script below sets every
  sample's size factor to 1) rather than normalizing twice. DESeq2 still
  performs its own dispersion estimation and negative-binomial significance
  testing on top of that.

  Reported overdispersion is the mean overdispersion across the
  NON-SELECTED replicates only (a simple (Var-Mean)/Mean^2 estimate,
  computed before DESeq2 runs) - not DESeq2's own per-group dispersion
  estimate, which this program does not report.

  NOTE: JAGS is no longer used anywhere in this pipeline. The original
  DEEPN's Bayesian JAGS-based statistics have been fully replaced by this
  DESeq2-based approach.

R SCRIPT (run_Y2H_enrichement_stats.R)
----------------------------------------------------------------------
"""

JUNCTION_COLUMNS = ['inframe_inorf', 'upstream', 'in_orf', 'downstream', 'in_frame', 'backwards', 'intron']
# get_junction_stats() keys its 'in-frame and in-ORF' percentage as 'frame_orf';
# everywhere else the column header matches the dict key exactly.
JUNCTION_COLUMN_KEYS = {'inframe_inorf': 'frame_orf'}


class MyListWidgetItem(QtWidgets.QListWidgetItem):
    def __init__(self, *args):
        super(MyListWidgetItem, self).__init__(*args)
        self.data = ''

    def GetData(self):
        return self.data


class RInstallThread(QtCore.QThread):
    lineOutput = QtCore.pyqtSignal(str)
    finished_ok = QtCore.pyqtSignal(bool, str)

    def __init__(self, rscript_exe, packages, parent=None):
        QtCore.QThread.__init__(self, parent)
        self.rscript_exe = rscript_exe
        self.packages = packages

    def run(self):
        try:
            sc.install_r_packages(self.rscript_exe, self.packages, progress_callback=self.lineOutput.emit)
            self.finished_ok.emit(True, '')
        except sc.RScriptError as e:
            self.finished_ok.emit(False, str(e))


class CollateThread(QtCore.QThread):
    lineOutput = QtCore.pyqtSignal(str)
    finished_ok = QtCore.pyqtSignal(bool, str)

    def __init__(self, parent_gui):
        QtCore.QThread.__init__(self, parent_gui)
        self.g = parent_gui

    def run(self):
        try:
            self.g._do_collate_and_stats(self.lineOutput.emit)
            self.finished_ok.emit(True, '')
        except Exception as e:
            self.finished_ok.emit(False, str(e))


class Stat_Maker_Gui(QtWidgets.QMainWindow, form_class):
    def __init__(self, *args):
        super(Stat_Maker_Gui, self).__init__(*args)
        self.setupUi(self)

        self.vec_sel_list.dropped.connect(self.vec_sel_fileDropped)
        self.vec_nonsel_list.dropped.connect(self.vec_nonsel_fileDropped)
        self.vec_sel_list_2.dropped.connect(self.vec_sel_fileDropped_2)
        self.vec_nonsel_list_2.dropped.connect(self.vec_nonsel_fileDropped_2)
        self.bait_sel_list.dropped.connect(self.bait_sel_fileDropped)
        self.bait_nonsel_list.dropped.connect(self.bait_nonsel_fileDropped)

        self.vec_sel_list.deleted.connect(self.file_deleted)
        self.vec_nonsel_list.deleted.connect(self.file_deleted)
        self.vec_sel_list_2.deleted.connect(self.file_deleted)
        self.vec_nonsel_list_2.deleted.connect(self.file_deleted)
        self.bait_sel_list.deleted.connect(self.file_deleted)
        self.bait_nonsel_list.deleted.connect(self.file_deleted)

        self.folder_choice_btn.clicked.connect(self.on_folder_choice_btn_clicked)
        self.verify_r_btn.clicked.connect(self.on_verify_r_btn_clicked)
        self.collate_btn.clicked.connect(self.on_collate_btn_clicked)
        self.method_info_btn.clicked.connect(self.on_method_info_btn_clicked)
        self.quit_btn.clicked.connect(self.on_quit_btn_clicked)

        self.fileio = f.fileio()
        self.data = {}
        self.directory = '~/'
        self.rscript_exe = sc.find_rscript()
        self.r_verified = False
        self.r_thread = None
        self.collate_thread = None

        self.log("Stat Maker v2 - collation + DESeq2 statistics")
        if self.rscript_exe:
            self.log("Found Rscript: %s" % self.rscript_exe)
        else:
            self.log("Rscript not found on this machine yet. Click 'Verify R Installation'.")

    def log(self, text):
        self.log_text.append(text)
        self.log_text.ensureCursorVisible()

    # ---------- R environment check/install ----------

    @QtCore.pyqtSlot()
    def on_verify_r_btn_clicked(self):
        self.rscript_exe = sc.find_rscript()
        if self.rscript_exe is None:
            box = QtWidgets.QMessageBox(self)
            box.setWindowTitle("R Not Found")
            box.setText("R does not appear to be installed on this Mac.\n\n"
                        "DEEPN's statistics step needs R (with DESeq2) to run.")
            download_btn = box.addButton("Download R...", QtWidgets.QMessageBox.ActionRole)
            box.addButton("Cancel", QtWidgets.QMessageBox.RejectRole)
            box.exec_()
            if box.clickedButton() == download_btn:
                QtGui.QDesktopServices.openUrl(QtCore.QUrl(sc.CRAN_DOWNLOAD_URL))
                self.log("Opened %s - after installing R, click 'Verify R Installation' again." % sc.CRAN_DOWNLOAD_URL)
            return

        self.log("Rscript found: %s" % self.rscript_exe)
        self.log("Checking required R packages...")
        self.verify_r_btn.setEnabled(False)
        QtWidgets.QApplication.processEvents()
        status = sc.check_r_packages(self.rscript_exe)
        missing = [p for p, ok in status.items() if not ok]
        for p, ok in status.items():
            self.log("  %-10s %s" % (p, "OK" if ok else "MISSING"))

        if not missing:
            self.r_verified = True
            self.log("All required R packages are installed.")
            self.verify_r_btn.setText("R Installation Verified")
            self.verify_r_btn.setEnabled(False)
            self._update_collate_enabled()
            return

        reply = QtWidgets.QMessageBox.question(
            self, "Install R Packages?",
            "The following R packages need to be installed:\n\n%s\n\n"
            "This can take several minutes (downloading/compiling DESeq2 and its dependencies). Proceed?"
            % ", ".join(missing),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No
        )
        if reply != QtWidgets.QMessageBox.Yes:
            self.verify_r_btn.setEnabled(True)
            return

        self.log("Installing: %s ..." % ", ".join(missing))
        self.r_thread = RInstallThread(self.rscript_exe, missing, self)
        self.r_thread.lineOutput.connect(self.log)
        self.r_thread.finished_ok.connect(self._r_install_finished)
        self.r_thread.start()

    @QtCore.pyqtSlot(bool, str)
    def _r_install_finished(self, ok, err):
        self.verify_r_btn.setEnabled(True)
        if not ok:
            self.log("R package installation FAILED: %s" % err)
            QtWidgets.QMessageBox.critical(self, "Installation Failed", err[:2000])
            return
        self.log("R package installation finished. Re-checking...")
        status = sc.check_r_packages(self.rscript_exe)
        missing = [p for p, ok2 in status.items() if not ok2]
        if missing:
            self.log("Still missing: %s" % ", ".join(missing))
        else:
            self.r_verified = True
            self.log("All required R packages are installed.")
            self.verify_r_btn.setText("R Installation Verified")
            self.verify_r_btn.setEnabled(False)
            self._update_collate_enabled()

    # ---------- folder / file selection ----------

    def initialize_folders(self, directory):
        self.fileio.create_new_folder(directory, "stat_maker_output")

    @QtCore.pyqtSlot()
    def on_folder_choice_btn_clicked(self):
        directory = str(QtWidgets.QFileDialog.getExistingDirectory(QtWidgets.QFileDialog(), "Locate Work Folder",
                                                               os.path.expanduser("~"),
                                                               QtWidgets.QFileDialog.ShowDirsOnly))
        if not directory:
            return
        self.directory = directory
        self.initialize_folders(self.directory)
        self.load_gene_summary_files()

    def load_gene_summary_files(self):
        try:
            dirlist = os.listdir(os.path.join(self.directory, 'gene_count_summary'))
            self.file_list.clear()
            for file in sorted(dirlist):
                if not re.match(r'^\.', file) and re.match(r'.+summary\.csv', file):
                    path = os.path.join(self.directory, 'gene_count_summary', file)
                    iconProvider = QtWidgets.QFileIconProvider()
                    fileInfo = QtCore.QFileInfo(path)
                    icon = iconProvider.icon(fileInfo)
                    item = MyListWidgetItem()
                    item.data = path
                    item.setIcon(icon)
                    item.setText(file)
                    self.file_list.addItem(item)
        except OSError:
            pass

    def fileDropped(self, list_widget, url, key):
        if os.path.exists(url):
            iconProvider = QtWidgets.QFileIconProvider()
            fileInfo = QtCore.QFileInfo(url)
            icon = iconProvider.icon(fileInfo)
            item = MyListWidgetItem()
            item.data = url
            item.setIcon(icon)
            item.setText(os.path.basename(os.path.normpath(url)))
            list_widget.clear()
            list_widget.addItem(item)
            self.data[key] = url
        self._update_collate_enabled()

    def file_deleted(self, path):
        for key in list(self.data.keys()):
            if self.data[key] == path:
                del self.data[key]
        self._update_collate_enabled()

    def vec_sel_fileDropped(self, path):
        self.fileDropped(self.vec_sel_list, str(path), 'Vector1S')

    def vec_nonsel_fileDropped(self, path):
        self.fileDropped(self.vec_nonsel_list, str(path), 'Vector1N')

    def vec_sel_fileDropped_2(self, path):
        self.fileDropped(self.vec_sel_list_2, str(path), 'Vector2S')

    def vec_nonsel_fileDropped_2(self, path):
        self.fileDropped(self.vec_nonsel_list_2, str(path), 'Vector2N')

    def bait_sel_fileDropped(self, path):
        self.fileDropped(self.bait_sel_list, str(path), 'BaitS')

    def bait_nonsel_fileDropped(self, path):
        self.fileDropped(self.bait_nonsel_list, str(path), 'BaitN')

    def _update_collate_enabled(self):
        required = ('Vector1S', 'Vector1N', 'Vector2S', 'Vector2N', 'BaitS', 'BaitN')
        ready = all(k in self.data for k in required)
        self.collate_btn.setEnabled(ready)

    # ---------- collate + stats ----------

    @QtCore.pyqtSlot()
    def on_collate_btn_clicked(self):
        if not self.r_verified:
            reply = QtWidgets.QMessageBox.question(
                self, "R Not Verified",
                "R installation hasn't been verified yet this session. Continue anyway "
                "(the statistics step will fail if R/DESeq2 isn't actually available)?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if reply != QtWidgets.QMessageBox.Yes:
                return

        self.collate_btn.setEnabled(False)
        self.statusbar.showMessage("Running collation + statistics...")
        self.collate_thread = CollateThread(self)
        self.collate_thread.lineOutput.connect(self.log)
        self.collate_thread.finished_ok.connect(self._collate_finished)
        self.collate_thread.start()

    @QtCore.pyqtSlot(bool, str)
    def _collate_finished(self, ok, err):
        self._update_collate_enabled()
        if ok:
            self.statusbar.showMessage("Done.", 8000)
        else:
            self.statusbar.showMessage("Failed - see log.", 8000)
            self.log("ERROR: %s" % err)
            QtWidgets.QMessageBox.critical(self, "Collation Failed", err[:2000])

    def _do_collate_and_stats(self, log_cb):
        directory = self.directory
        out_folder = os.path.join(directory, 'stat_maker_output')
        self.fileio.create_new_folder(directory, 'stat_maker_output')

        csv_paths = {k: self.data[k] for k in ('Vector1S', 'Vector1N', 'Vector2S', 'Vector2N', 'BaitS', 'BaitN')}

        log_cb("Checking gene counts are consistent across all 6 files...")
        consistent, counts = sc.check_gene_counts_consistent(list(csv_paths.values()))
        for k, path in csv_paths.items():
            log_cb("  %-10s %6d genes  (%s)" % (k, counts[path], os.path.basename(path)))
        if not consistent:
            log_cb("WARNING: gene counts differ across files - they may not have been "
                   "processed with the same gene dictionary/list.")

        log_cb("Loading PPM data (summed across scaffold-duplicate gene rows)...")
        ppm_data = {k: sc.load_ppm_summed(path) for k, path in csv_paths.items()}

        timestamp = datetime.now().strftime('%d%b%Y-%H%M%S')
        run_out_dir = os.path.join(out_folder, 'run_%s' % timestamp)
        os.makedirs(run_out_dir, exist_ok=True)
        collated_path = os.path.join(run_out_dir, 'collated_input.csv')
        sc.build_collated_input(ppm_data, run_out_dir, collated_path)
        log_cb("Wrote %s" % collated_path)

        log_cb("Running DESeq2 statistics (this can take a minute)...")
        r_output, overdispersion = sc.run_r_script(R_SCRIPT_PATH, collated_path, self.rscript_exe)
        for line in r_output.splitlines():
            log_cb("  [R] " + line)
        if overdispersion is not None:
            log_cb("Non-selected overdispersion: %.4f" % overdispersion)

        stats_csv = os.path.join(run_out_dir, 'everything_combined.csv')
        stats = sc.load_stats_output(stats_csv)

        log_cb("Loading .bqp junction data for collation...")
        bqp_data = {}
        for k, path in csv_paths.items():
            bqp_name = os.path.basename(path).replace('_summary.csv', '') + '.bqp'
            bqp_path = os.path.join(directory, 'blast_results_query', bqp_name)
            if os.path.exists(bqp_path):
                bqp_data[k] = sc.load_bqp(bqp_path)
            else:
                log_cb("  WARNING: no .bqp found for %s (expected %s) - junction stats will be 0 for it" % (k, bqp_name))
                bqp_data[k] = {}

        dataset_order = ['BaitN', 'BaitS', 'Vector1N', 'Vector2N', 'Vector1S', 'Vector2S']
        dataset_labels = {
            'BaitN': 'Bait_Non-Selected', 'BaitS': 'Bait_Selected',
            'Vector1N': 'Vector_Non-Selected_1', 'Vector2N': 'Vector_Non-Selected_2',
            'Vector1S': 'Vector_Selected_1', 'Vector2S': 'Vector_Selected_2',
        }

        genes = sorted(ppm_data['BaitS'].keys())
        log_cb("Collating junction stats for %d genes..." % len(genes))

        out_name = 'stat_maker_%s.xlsx' % timestamp
        out_path = os.path.join(out_folder, out_name)
        workbook = xls.Workbook(out_path)
        ws = workbook.add_worksheet('results')

        row = 0
        for k in ('BaitS', 'BaitN', 'Vector1S', 'Vector1N', 'Vector2S', 'Vector2N'):
            ws.write(row, 0, dataset_labels[k])
            ws.write(row, 1, csv_paths[k])
            row += 2
        if overdispersion is not None:
            ws.write(row, 0, 'Overdispersion (non-selected)')
            ws.write(row, 1, overdispersion)
            row += 2

        header_row = row
        headers = ['Gene', 'pvalue_bait', 'log2FoldChange_bait', 'Enrichment_score_bait',
                  'pvalue_vector', 'log2FoldChange_vector', 'Enrichment_score_vector']
        for c, h in enumerate(headers):
            ws.write(header_row + 1, c, h)
        block_col = len(headers)
        for k in dataset_order:
            ws.write(header_row, block_col, dataset_labels[k])
            for j, col_name in enumerate(JUNCTION_COLUMNS):
                ws.write(header_row + 1, block_col + j, col_name)
            block_col += len(JUNCTION_COLUMNS)

        data_row = header_row + 2
        for gene in genes:
            col = 0
            ws.write(data_row, col, gene); col += 1
            srow = stats.get(gene, {})
            ws.write(data_row, col, srow.get('pvalue_bait', '')); col += 1
            ws.write(data_row, col, srow.get('log2FoldChange_bait', '')); col += 1
            ws.write(data_row, col, srow.get('Enrichment_score_bait', '')); col += 1
            ws.write(data_row, col, srow.get('pvalue_vector', '')); col += 1
            ws.write(data_row, col, srow.get('log2FoldChange_vector', '')); col += 1
            ws.write(data_row, col, srow.get('Enrichment_score_vector', '')); col += 1

            for k in dataset_order:
                gstat = sc.get_junction_stats(bqp_data[k], gene)
                for jc in JUNCTION_COLUMNS:
                    ws.write(data_row, col, gstat[JUNCTION_COLUMN_KEYS.get(jc, jc)])
                    col += 1
            data_row += 1

        workbook.close()
        log_cb("Wrote %s" % out_path)

    @QtCore.pyqtSlot()
    def on_method_info_btn_clicked(self):
        try:
            with open(R_SCRIPT_PATH, encoding='utf-8') as fh:
                r_source = fh.read()
        except (OSError, UnicodeDecodeError) as e:
            r_source = "(Could not read %s: %s)" % (R_SCRIPT_PATH, e)

        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("Method Info")
        dialog.resize(800, 700)
        layout = QtWidgets.QVBoxLayout(dialog)
        text = QtWidgets.QPlainTextEdit(dialog)
        text.setReadOnly(True)
        text.setFont(QtGui.QFont("Menlo", 12))
        text.setPlainText(METHOD_INFO_TEXT + r_source)
        layout.addWidget(text)
        close_btn = QtWidgets.QPushButton("Close", dialog)
        close_btn.clicked.connect(dialog.accept)
        layout.addWidget(close_btn)
        dialog.exec_()

    @QtCore.pyqtSlot()
    def on_quit_btn_clicked(self):
        app.quit()
        sys.exit()


def appExit():
    app.quit()
    sys.exit()

if __name__ == '__main__':
    form = Stat_Maker_Gui()
    form.show()
    app.aboutToQuit.connect(appExit)
    app.exec_()
