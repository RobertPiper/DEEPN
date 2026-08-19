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
form_class, base_class = uic.loadUiType(os.path.join('ui', 'Stat_Maker_v3.ui'))

R_SCRIPT_PATH = os.path.join('r_scripts', 'run_Y2H_enrichement_stats_v3.R')
BAYESIAN_R_SCRIPT_PATH = os.path.join('r_scripts', 'run_Y2H_bayesian_stats_v3.R')
BAYESIAN_THRESHOLD_PPM = 3  # matches Stat Maker v1/v2's threshold_sbx default

METHOD_INFO_TEXT = """STAT MAKER v3 - METHOD OVERVIEW

WHAT THIS PROGRAM DOES

Stat Maker compares one or two Y2H baits (each Selected vs Non-Selected)
against a Vector control (2 replicates, Selected vs Non-Selected each),
then combines the results into one table, one row per gene:

  1. STATISTICS (via DESeq2, see R code below)
     For each gene: pvalue, log2FoldChange, and a rank-based Enrichment_score,
     computed separately for Bait1 (Selected vs Non-Selected), Vector
     (Selected vs Non-Selected), and - if a Bait2 dataset is also
     supplied - for Bait2 as well, plus a Bait1-vs-Bait2 specificity
     contrast (pvalue_specificity/log2FoldChange_specificity: positive
     means more enriched under Bait1 selection than Bait2 selection, i.e.
     more specific to Bait1; negative means the opposite).

     Bait2 is optional. With just Bait1 + Vector supplied, this produces
     the exact same Bait1/Vector columns Stat Maker v2 always did (just
     named "_bait1" instead of "_bait") - no specificity columns, since
     there's no second bait to compare against.

  2. JUNCTION COLLATION (from BLAST results, .bqp files)
     For each gene, in each input dataset, the percentage of junctions
     that are:
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

  Adapted from the method described in:

    Velasquez-Zapata V, Elmore JM, Banerjee S, Dorman KS, Wise RP (2021)
    Next-generation yeast-two-hybrid analysis with Y2H-SCORES identifies
    novel interactors of the MLA immune receptor. PLOS Computational
    Biology 17(4): e1008890. https://doi.org/10.1371/journal.pcbi.1008890

  The enrichment contrast (Selected vs Non-Selected per bait/vector) is
  the same streamlined single-file approach Stat Maker v2 used. The
  specificity contrast (Bait1-selected vs Bait2-selected, when Bait2 is
  supplied) is a direct pairwise DESeq2 comparison - a scaled-down version
  of Y2H-SCORES' own specificity score (Software/spec_score.R), which
  generalizes to any number of baits via pairwise combinations; this
  build only ever has two baits, so there's only ever one pair to compare,
  and none of that generalization (or its kde2d-based fold-change
  reweighting) is needed here.

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

STATISTICS METHOD (Bayesian/MCMC, via JAGS)

  This is DEEPN's original statistics method (Patrick Breheny,
  pbreheny/deepn), run here from an unmodified, vendored copy of that
  package (r_scripts/deepn_bayesian_src/) rather than re-derived. It
  requires JAGS 4.x installed separately (see "Verify R + JAGS
  Installation") - CRAN's rjags package hard-rejects JAGS 5.x.

  A per-gene negative-binomial hierarchical model is fit by MCMC
  (4 chains, 1500 adapt/burnin/sample iterations each - runMCMC() in the
  vendored mcmc.R). AdjEnr1/AdjEnr2 are posterior median log2 effect
  sizes (Bait1/Bait2 vs Vector) and reproduce almost exactly run-to-run.
  The p-columns (pBait1_Vec, pBait2_Vec, pBait1_Bait2, pBait2_Bait1,
  pBait1, pBait2) are pnorm(median/MAD) - a posterior probability, not a
  frequentist p-value - and because pnorm saturates hard toward 0/1, these
  columns can flip between runs for genes whose effect is only moderately
  confident, even though the underlying AdjEnr estimate barely moves. This
  is a known, inherent property of the original method (not something
  introduced in this port) - treat the p-columns as directional/qualitative
  and AdjEnr as the more stable quantity when the two seem to disagree.

  With just Bait1 + Vector supplied (no Bait2), the single-bait model
  (sModel1.jag) is used instead, reporting just AdjEnr/p.

R SCRIPT (run_Y2H_enrichement_stats_v3.R)
----------------------------------------------------------------------
"""

JUNCTION_COLUMNS = ['inframe_inorf', 'upstream', 'in_orf', 'downstream', 'in_frame', 'backwards', 'intron']
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


class JagsDownloadThread(QtCore.QThread):
    lineOutput = QtCore.pyqtSignal(str)
    finished_ok = QtCore.pyqtSignal(bool, str)

    def run(self):
        try:
            pkg_path = sc.download_jags_pkg(progress_callback=self.lineOutput.emit)
            self.finished_ok.emit(True, pkg_path)
        except sc.RScriptError as e:
            self.finished_ok.emit(False, str(e))


class DeepnInstallThread(QtCore.QThread):
    lineOutput = QtCore.pyqtSignal(str)
    finished_ok = QtCore.pyqtSignal(bool, str)

    def __init__(self, rscript_exe, parent=None):
        QtCore.QThread.__init__(self, parent)
        self.rscript_exe = rscript_exe

    def run(self):
        try:
            sc.install_deepn_package(self.rscript_exe, progress_callback=self.lineOutput.emit)
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
        self.bait1_sel_list.dropped.connect(self.bait1_sel_fileDropped)
        self.bait1_nonsel_list.dropped.connect(self.bait1_nonsel_fileDropped)
        self.bait2_sel_list.dropped.connect(self.bait2_sel_fileDropped)
        self.bait2_nonsel_list.dropped.connect(self.bait2_nonsel_fileDropped)

        self.vec_sel_list.deleted.connect(self.file_deleted)
        self.vec_nonsel_list.deleted.connect(self.file_deleted)
        self.vec_sel_list_2.deleted.connect(self.file_deleted)
        self.vec_nonsel_list_2.deleted.connect(self.file_deleted)
        self.bait1_sel_list.deleted.connect(self.file_deleted)
        self.bait1_nonsel_list.deleted.connect(self.file_deleted)
        self.bait2_sel_list.deleted.connect(self.file_deleted)
        self.bait2_nonsel_list.deleted.connect(self.file_deleted)

        # Handler names deliberately do NOT match the on_<widget>_<signal>
        # pattern - setupUi()'s connectSlotsByName() binds a matching name to
        # BOTH of clicked's signal overloads (clicked() and clicked(bool)),
        # firing the handler twice per click even with no explicit connect
        # at all. Avoiding the naming convention and connecting explicitly,
        # once, is the only reliable way to get exactly one call per click.
        self.folder_choice_btn.clicked.connect(self.folder_choice_clicked)
        self.verify_env_btn.clicked.connect(self.verify_env_clicked)
        self.collate_btn.clicked.connect(self.collate_clicked)
        self.method_info_btn.clicked.connect(self.method_info_clicked)
        self.quit_btn.clicked.connect(self.quit_clicked)

        self.fileio = f.fileio()
        self.data = {}
        self.directory = '~/'
        self.rscript_exe = sc.find_rscript()
        self.jags_exe = sc.find_jags()
        self.jags_compatible = False
        if self.jags_exe:
            self.jags_compatible, self.jags_version = sc.check_jags_version(self.jags_exe)
        else:
            self.jags_version = None
        self.env_verified = False
        self.r_thread = None
        self.jags_download_thread = None
        self.deepn_thread = None
        self.collate_thread = None

        self.log("Stat Maker v3 - collation + DESeq2 statistics (Bait1 required, Bait2 optional)")
        self.log("Found Rscript: %s" % self.rscript_exe if self.rscript_exe else
                 "Rscript not found on this machine yet.")
        if self.jags_exe and self.jags_compatible:
            self.log("Found jags %s (compatible): %s" % (self.jags_version, self.jags_exe))
        elif self.jags_exe:
            self.log("Found jags %s at %s, but rjags requires JAGS 4.x - this version won't work." %
                     (self.jags_version, self.jags_exe))
        else:
            self.log("jags not found on this machine yet.")
        if not (self.rscript_exe and self.jags_exe and self.jags_compatible):
            self.log("Click 'Verify R + JAGS Installation'.")

    def log(self, text):
        self.log_text.append(text)
        self.log_text.ensureCursorVisible()

    # ---------- R + JAGS environment check/install ----------
    # Three tiers, checked in order: R itself and JAGS itself are separate
    # standalone programs neither can install for the other - if missing,
    # this only ever points the user at a manual install (CRAN's own
    # installer for R; the official SourceForge installer for JAGS - NOT
    # Homebrew, whose `jags` formula is 5.0.0, a version CRAN's rjags
    # package hard-rejects at its own configure step). Once both exist,
    # the R *packages* each needs (DESeq2 etc., and now rjags/runjags) are
    # auto-installable through R's own package manager, same as always.

    @QtCore.pyqtSlot()
    def verify_env_clicked(self):
        self.verify_env_btn.setEnabled(False)
        self._verify_rscript_step()

    def _verify_rscript_step(self):
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
                self.log("Opened %s - after installing R, click 'Verify R + JAGS Installation' again." % sc.CRAN_DOWNLOAD_URL)
            self.verify_env_btn.setEnabled(True)
            return

        self.log("Rscript found: %s" % self.rscript_exe)
        self._verify_jags_step()

    def _verify_jags_step(self):
        self.jags_exe = sc.find_jags()
        self.jags_compatible = False
        self.jags_version = None
        if self.jags_exe:
            self.jags_compatible, self.jags_version = sc.check_jags_version(self.jags_exe)

        if not self.jags_compatible:
            box = QtWidgets.QMessageBox(self)
            box.setWindowTitle("JAGS 4.x Required")
            if self.jags_exe:
                box.setText(
                    "Found JAGS %s at %s, but the Bayesian statistics method needs "
                    "JAGS 4.x specifically - rjags (the R/JAGS bridge) rejects JAGS 5.x.\n\n"
                    "Note: Homebrew's 'jags' formula currently only builds 5.x and is NOT "
                    "usable here. Install the official JAGS 4.3.2 package below instead."
                    % (self.jags_version, self.jags_exe))
            else:
                box.setText(
                    "JAGS (the MCMC engine the Bayesian statistics method needs) does not "
                    "appear to be installed on this Mac.\n\n"
                    "Install JAGS 4.3.2 from the official installer below - NOT Homebrew's "
                    "'jags' formula, which is version 5.x and incompatible with rjags.")
            fetch_btn = box.addButton("Download and Open Installer", QtWidgets.QMessageBox.ActionRole)
            download_btn = box.addButton("Open Download Page", QtWidgets.QMessageBox.ActionRole)
            box.addButton("Cancel", QtWidgets.QMessageBox.RejectRole)
            box.exec_()
            clicked = box.clickedButton()
            if clicked == fetch_btn:
                self.log("Downloading %s ..." % sc.JAGS_PKG_URL)
                self.jags_download_thread = JagsDownloadThread(self)
                self.jags_download_thread.lineOutput.connect(self.log)
                self.jags_download_thread.finished_ok.connect(self._jags_download_finished)
                self.jags_download_thread.start()
                return
            elif clicked == download_btn:
                QtGui.QDesktopServices.openUrl(QtCore.QUrl(sc.JAGS_DOWNLOAD_URL))
                self.log("Opened %s - after installing JAGS 4.x, click 'Verify R + JAGS Installation' again." % sc.JAGS_DOWNLOAD_URL)
            self.verify_env_btn.setEnabled(True)
            return

        self.log("jags %s found (compatible): %s" % (self.jags_version, self.jags_exe))
        self._verify_packages_step()

    @QtCore.pyqtSlot(bool, str)
    def _jags_download_finished(self, ok, result):
        if not ok:
            self.log("JAGS download FAILED: %s" % result)
            QtWidgets.QMessageBox.critical(self, "Download Failed", result[:2000])
            self.verify_env_btn.setEnabled(True)
            return
        pkg_path = result
        self.log("Opening %s in Installer.app - follow the prompts and enter your Mac "
                 "password when asked, then click 'Verify R + JAGS Installation' again." % pkg_path)
        sc.open_installer(pkg_path)
        self.verify_env_btn.setEnabled(True)

    def _verify_packages_step(self):
        self.log("Checking required R packages...")
        QtWidgets.QApplication.processEvents()
        packages = sc.REQUIRED_R_PACKAGES + sc.BAYESIAN_R_PACKAGES
        status = sc.check_r_packages(self.rscript_exe, packages)
        missing = [p for p, ok in status.items() if not ok]
        for p, ok in status.items():
            self.log("  %-10s %s" % (p, "OK" if ok else "MISSING"))

        if not missing:
            self._verify_deepn_step()
            return

        reply = QtWidgets.QMessageBox.question(
            self, "Install R Packages?",
            "The following R packages need to be installed:\n\n%s\n\n"
            "This can take several minutes (downloading/compiling DESeq2 and its dependencies). Proceed?"
            % ", ".join(missing),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No
        )
        if reply != QtWidgets.QMessageBox.Yes:
            self.verify_env_btn.setEnabled(True)
            return

        self.log("Installing: %s ..." % ", ".join(missing))
        self.r_thread = RInstallThread(self.rscript_exe, missing, self)
        self.r_thread.lineOutput.connect(self.log)
        self.r_thread.finished_ok.connect(self._r_install_finished)
        self.r_thread.start()

    @QtCore.pyqtSlot(bool, str)
    def _r_install_finished(self, ok, err):
        if not ok:
            self.log("R package installation FAILED: %s" % err)
            QtWidgets.QMessageBox.critical(self, "Installation Failed", err[:2000])
            self.verify_env_btn.setEnabled(True)
            return
        self.log("R package installation finished. Re-checking...")
        packages = sc.REQUIRED_R_PACKAGES + sc.BAYESIAN_R_PACKAGES
        status = sc.check_r_packages(self.rscript_exe, packages)
        missing = [p for p, ok2 in status.items() if not ok2]
        if missing:
            self.log("Still missing: %s" % ", ".join(missing))
            self.verify_env_btn.setEnabled(True)
        else:
            self._verify_deepn_step()

    def _verify_deepn_step(self):
        self.log("Checking vendored 'deepn' package (pbreheny/deepn, hardwired local copy)...")
        QtWidgets.QApplication.processEvents()
        if sc.check_deepn_package(self.rscript_exe):
            self.log("  deepn      OK")
            self._env_fully_verified()
            return

        self.log("  deepn      MISSING - installing from vendored local source...")
        self.deepn_thread = DeepnInstallThread(self.rscript_exe, self)
        self.deepn_thread.lineOutput.connect(self.log)
        self.deepn_thread.finished_ok.connect(self._deepn_install_finished)
        self.deepn_thread.start()

    @QtCore.pyqtSlot(bool, str)
    def _deepn_install_finished(self, ok, err):
        if not ok:
            self.log("deepn package installation FAILED: %s" % err)
            QtWidgets.QMessageBox.critical(self, "Installation Failed", err[:2000])
            self.verify_env_btn.setEnabled(True)
            return
        self.log("deepn package installed.")
        self._env_fully_verified()

    def _env_fully_verified(self):
        self.env_verified = True
        self.log("R + JAGS environment fully verified.")
        self.verify_env_btn.setText("R + JAGS Verified")
        self.verify_env_btn.setEnabled(False)
        self._update_collate_enabled()

    # ---------- folder / file selection ----------

    def initialize_folders(self, directory):
        self.fileio.create_new_folder(directory, "stat_maker_output")

    @QtCore.pyqtSlot()
    def folder_choice_clicked(self):
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
        else:
            self.data.pop(key, None)
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

    def bait1_sel_fileDropped(self, path):
        self.fileDropped(self.bait1_sel_list, str(path), 'Bait1S')

    def bait1_nonsel_fileDropped(self, path):
        self.fileDropped(self.bait1_nonsel_list, str(path), 'Bait1N')

    def bait2_sel_fileDropped(self, path):
        self.fileDropped(self.bait2_sel_list, str(path), 'Bait2S')

    def bait2_nonsel_fileDropped(self, path):
        self.fileDropped(self.bait2_nonsel_list, str(path), 'Bait2N')

    def _has_bait2(self):
        return 'Bait2S' in self.data and 'Bait2N' in self.data

    def _update_collate_enabled(self):
        required = ('Vector1S', 'Vector1N', 'Vector2S', 'Vector2N', 'Bait1S', 'Bait1N')
        ready = all(k in self.data for k in required)
        # Bait2 is optional, but if only one of its two files is dropped
        # (not both), that's a half-finished state - don't run until it's
        # either both-present or both-absent.
        bait2_keys = ('Bait2S', 'Bait2N')
        bait2_partial = any(k in self.data for k in bait2_keys) and not all(k in self.data for k in bait2_keys)
        self.collate_btn.setEnabled(ready and not bait2_partial)

    # ---------- collate + stats ----------

    @QtCore.pyqtSlot()
    def collate_clicked(self):
        if not self.env_verified:
            reply = QtWidgets.QMessageBox.question(
                self, "Environment Not Verified",
                "R + JAGS installation hasn't been verified yet this session. Continue anyway "
                "(the statistics step will fail if R/DESeq2/JAGS isn't actually available)?",
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

        has_bait2 = self._has_bait2()
        keys = ['Vector1S', 'Vector1N', 'Vector2S', 'Vector2N', 'Bait1S', 'Bait1N']
        if has_bait2:
            keys += ['Bait2S', 'Bait2N']
        csv_paths = {k: self.data[k] for k in keys}

        log_cb("Checking gene counts are consistent across all %d files..." % len(csv_paths))
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
        sc.build_collated_input_v3(ppm_data, run_out_dir, collated_path, has_bait2)
        log_cb("Wrote %s%s" % (collated_path, " (Bait1 + Bait2 + Vector)" if has_bait2 else " (Bait1 + Vector only)"))

        log_cb("Running DESeq2 statistics (this can take a minute)...")
        r_output, overdispersion = sc.run_r_script(R_SCRIPT_PATH, collated_path, self.rscript_exe)
        for line in r_output.splitlines():
            log_cb("  [R] " + line)
        if overdispersion is not None:
            log_cb("Non-selected overdispersion: %.4f" % overdispersion)

        stats_csv = os.path.join(run_out_dir, 'everything_combined.csv')
        stats = sc.load_stats_output(stats_csv)

        log_cb("Running Bayesian/MCMC statistics (analyzeDeepn via JAGS - this is the slow "
               "step, several minutes)...")
        try:
            bayes_output, bayes_csv, bayes_overdispersion = sc.run_bayesian_r_script(
                BAYESIAN_R_SCRIPT_PATH, csv_paths, BAYESIAN_THRESHOLD_PPM, run_out_dir, self.rscript_exe)
            for line in bayes_output.splitlines():
                log_cb("  [R/JAGS] " + line)
            bayesian_stats = sc.load_bayesian_stats_output(bayes_csv)
        except sc.RScriptError as e:
            log_cb("Bayesian/MCMC statistics FAILED (DESeq2 results above are still valid): %s" % e)
            bayesian_stats = {}
            bayes_overdispersion = {}

        log_cb("Loading .bqp junction data for collation...")
        bqp_data = {}
        for k, path in csv_paths.items():
            basename = os.path.basename(path).replace('_summary.csv', '')
            bqp_path = sc.find_bqp_path(directory, basename)
            if bqp_path:
                if os.path.basename(bqp_path) == basename + '.bqp':
                    log_cb("  %-10s using legacy-named 5p file (no prefix): %s" % (k, os.path.basename(bqp_path)))
                bqp_data[k] = sc.load_bqp(bqp_path)
            else:
                log_cb("  WARNING: no .bqp found for %s (checked 5p_%s.bqp and %s.bqp) - junction stats will be 0 for it"
                       % (k, basename, basename))
                bqp_data[k] = {}

        dataset_order = ['Bait1N', 'Bait1S'] + (['Bait2N', 'Bait2S'] if has_bait2 else []) + \
                        ['Vector1N', 'Vector2N', 'Vector1S', 'Vector2S']
        dataset_labels = {
            'Bait1N': 'Bait1_Non-Selected', 'Bait1S': 'Bait1_Selected',
            'Bait2N': 'Bait2_Non-Selected', 'Bait2S': 'Bait2_Selected',
            'Vector1N': 'Vector_Non-Selected_1', 'Vector2N': 'Vector_Non-Selected_2',
            'Vector1S': 'Vector_Selected_1', 'Vector2S': 'Vector_Selected_2',
        }

        genes = sorted(ppm_data['Bait1S'].keys())
        log_cb("Collating junction stats for %d genes..." % len(genes))

        out_name = 'stat_maker_%s.xlsx' % timestamp
        out_path = os.path.join(out_folder, out_name)
        self._write_workbook(out_path, csv_paths, dataset_order, dataset_labels, overdispersion,
                             bayes_overdispersion, stats, bayesian_stats, bqp_data, genes, has_bait2, log_cb)
        log_cb("Wrote %s" % out_path)

    # ---------- output workbook ----------
    # Two-row column headers: an upper "group" row (merged, color-coded by
    # method - DESeq2 now, matching colors reserved for an MCMC/Bayesian
    # block once that method is added) sitting above the per-field header
    # row that the Excel autofilter buttons (sort/filter dropdowns, like a
    # pivot table) actually attach to.

    def _write_workbook(self, out_path, csv_paths, dataset_order, dataset_labels, overdispersion,
                        bayes_overdispersion, stats, bayesian_stats, bqp_data, genes, has_bait2, log_cb):
        workbook = xls.Workbook(out_path)
        ws = workbook.add_worksheet('results')

        fmt_group_deseq2 = workbook.add_format({'bold': True, 'align': 'center', 'valign': 'vcenter',
                                                'bg_color': '#D9E8FB', 'border': 1})
        fmt_group_bayesian = workbook.add_format({'bold': True, 'align': 'center', 'valign': 'vcenter',
                                                   'bg_color': '#E6D9F7', 'border': 1})
        fmt_group_junction = workbook.add_format({'bold': True, 'align': 'center', 'valign': 'vcenter',
                                                   'bg_color': '#E9E9E9', 'border': 1})
        fmt_field = workbook.add_format({'bold': True, 'border': 1})
        fmt_gene = workbook.add_format({'bold': True})
        fmt_section = workbook.add_format({'bold': True})

        row = 0
        for k in csv_paths:
            ws.write(row, 0, dataset_labels[k])
            ws.write(row, 1, csv_paths[k])
            row += 1

        # Both methods' overdispersion estimates, side by side for direct
        # comparison: analyzeDeepn()'s edgeR-based estimates (the original
        # Bayesian-method math, computed on raw counts) alongside DESeq2's
        # own estimate (computed on PPM) - deliberately kept separate since
        # they're different quantities from different pipelines, not
        # expected to match each other.
        overdisp_rows = [
            ('Vector non-select', bayes_overdispersion.get('vector_nonselect')),
            ('vector+Bait non-select', bayes_overdispersion.get('vector_bait_nonselect')),
            ('Vector select', bayes_overdispersion.get('vector_select')),
            ('DESeq Vector non-select', overdispersion),
        ]
        if any(v is not None for _, v in overdisp_rows):
            ws.write(row, 0, 'OVERDISPERSION', fmt_section)
            row += 1
            for label, value in overdisp_rows:
                ws.write(row, 0, label)
                if value is not None:
                    ws.write_number(row, 1, value)
                row += 1
        row += 1

        group_row = row
        field_row = row + 1

        # Gene column spans both header rows
        ws.merge_range(group_row, 0, field_row, 0, 'Gene', fmt_gene)

        col = 1
        deseq2_blocks = [('Bait1 Enrichment (DESeq2)', 'bait1', ['pvalue', 'log2FoldChange', 'Enrichment_score'])]
        if has_bait2:
            deseq2_blocks.append(('Bait2 Enrichment (DESeq2)', 'bait2', ['pvalue', 'log2FoldChange', 'Enrichment_score']))
        deseq2_blocks.append(('Vector Enrichment (DESeq2)', 'vector', ['pvalue', 'log2FoldChange', 'Enrichment_score']))
        if has_bait2:
            deseq2_blocks.append(('Bait1 vs Bait2 Specificity (DESeq2)', 'specificity', ['pvalue', 'log2FoldChange']))

        stat_field_cols = []  # (stats_dict_key, worksheet_col)
        for label, suffix, fields in deseq2_blocks:
            start_col = col
            for fld in fields:
                stats_key = '%s_%s' % (fld, suffix)
                ws.write(field_row, col, fld, fmt_field)
                stat_field_cols.append((stats_key, col))
                col += 1
            if col - 1 > start_col:
                ws.merge_range(group_row, start_col, group_row, col - 1, label, fmt_group_deseq2)
            else:
                ws.write(group_row, start_col, label, fmt_group_deseq2)

        # Bayesian/MCMC block (analyzeDeepn via JAGS). Column set mirrors
        # summary-psm.R's two branches: the 2-bait model reports AdjEnr1/
        # AdjEnr2 (posterior effect size per bait) plus four pairwise
        # probabilities and the two specificity-adjusted pBait1/pBait2; the
        # 1-bait model (no Bait2 supplied) reports just AdjEnr/p.
        if has_bait2:
            bayesian_fields = ['AdjEnr1', 'AdjEnr2', 'pBait1_Vec', 'pBait2_Vec',
                               'pBait1_Bait2', 'pBait2_Bait1', 'pBait1', 'pBait2']
        else:
            bayesian_fields = ['AdjEnr', 'p']
        bayesian_field_cols = []  # (bayesian_dict_key, worksheet_col)
        start_col = col
        for fld in bayesian_fields:
            ws.write(field_row, col, fld, fmt_field)
            bayesian_field_cols.append((fld, col))
            col += 1
        ws.merge_range(group_row, start_col, group_row, col - 1, 'Bayesian Enrichment (MCMC/JAGS)', fmt_group_bayesian)

        junction_field_cols = []  # (dataset_key, junction_col, worksheet_col)
        for k in dataset_order:
            start_col = col
            for jc in JUNCTION_COLUMNS:
                ws.write(field_row, col, jc, fmt_field)
                junction_field_cols.append((k, jc, col))
                col += 1
            ws.merge_range(group_row, start_col, group_row, col - 1, dataset_labels[k], fmt_group_junction)

        last_col = col - 1
        data_row = field_row + 1
        for gene in genes:
            ws.write(data_row, 0, gene)
            srow = stats.get(gene, {})
            for stats_key, wcol in stat_field_cols:
                # Leave the cell genuinely blank (not '') when a gene has no
                # value - Excel sorts numbers vs. text in reversed relative
                # order between ascending/descending, so a text '' placeholder
                # jumps to the top on a descending sort; a true blank cell
                # sorts to the end regardless of direction. Values also come
                # from csv.DictReader as strings - write as float so Excel
                # sorts/filters them numerically instead of lexicographically.
                # DESeq2 legitimately reports 'NA' for some genes (outlier/
                # independent filtering) - treat that the same as missing.
                if stats_key in srow and srow[stats_key] != 'NA':
                    ws.write_number(data_row, wcol, float(srow[stats_key]))
            brow = bayesian_stats.get(gene, {})
            for fld, wcol in bayesian_field_cols:
                if fld in brow and brow[fld] != 'NA':
                    ws.write_number(data_row, wcol, float(brow[fld]))
            for k, jc, wcol in junction_field_cols:
                gstat = sc.get_junction_stats(bqp_data[k], gene)
                ws.write(data_row, wcol, gstat[JUNCTION_COLUMN_KEYS.get(jc, jc)])
            data_row += 1

        # Sort/filter dropdowns on the field-name row, pivot-table style.
        ws.autofilter(field_row, 0, data_row - 1, last_col)
        ws.freeze_panes(data_row - len(genes), 1)
        workbook.close()

    @QtCore.pyqtSlot()
    def method_info_clicked(self):
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
    def quit_clicked(self):
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
