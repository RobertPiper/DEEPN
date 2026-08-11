#!/usr/bin/env python3
import os
import re
import sys
import threading
import glob
import pickle
import time
import json
import subprocess
import io
import traceback

from datetime import datetime
import functions.fileio_gui as f
import functions.stat_graph as statgraph
import functions.structures as strt
from PyQt5 import QtCore, QtGui, QtWidgets, uic
import DragDropListView

app = QtWidgets.QApplication(sys.argv)
form_class, base_class = uic.loadUiType(os.path.join('ui', 'Stat_Maker.ui'))


class MyListWidgetItem(QtWidgets.QListWidgetItem):
    def __init__(self, *args):
        super(MyListWidgetItem, self).__init__(*args)
        self.data = ''

    def GetData(self):
        return self.data


class Stat_Maker_Gui(QtWidgets.QMainWindow, form_class):
    def __init__(self, *args):
        super(Stat_Maker_Gui, self).__init__(*args)
        self.setupUi(self)
        self.vec_sel_list.dropped.connect(self.vec_sel_fileDropped)
        self.vec_nonsel_list.dropped.connect(self.vec_nonsel_fileDropped)
        self.bait1_sel_list.dropped.connect(self.bait1_sel_fileDropped)
        self.bait1_nonsel_list.dropped.connect(self.bait1_nonsel_fileDropped)
        self.bait2_sel_list.dropped.connect(self.bait2_sel_fileDropped)
        self.bait2_nonsel_list.dropped.connect(self.bait2_nonsel_fileDropped)

        self.vec_sel_list_2.dropped.connect(self.vec_sel_fileDropped_2)
        self.vec_nonsel_list_2.dropped.connect(self.vec_nonsel_fileDropped_2)
        self.bait1_sel_list_2.dropped.connect(self.bait1_sel_fileDropped_2)
        self.bait1_nonsel_list_2.dropped.connect(self.bait1_nonsel_fileDropped_2)
        self.bait2_sel_list_2.dropped.connect(self.bait2_sel_fileDropped_2)
        self.bait2_nonsel_list_2.dropped.connect(self.bait2_nonsel_fileDropped_2)

        self.vec_sel_list.deleted.connect(self.file_deleted)
        self.vec_nonsel_list.deleted.connect(self.file_deleted)
        self.bait1_sel_list.deleted.connect(self.file_deleted)
        self.bait1_nonsel_list.deleted.connect(self.file_deleted)
        self.bait2_sel_list.deleted.connect(self.file_deleted)
        self.bait2_nonsel_list.deleted.connect(self.file_deleted)

        self.vec_sel_list_2.deleted.connect(self.file_deleted)
        self.vec_nonsel_list_2.deleted.connect(self.file_deleted)
        self.bait1_sel_list_2.deleted.connect(self.file_deleted)
        self.bait1_nonsel_list_2.deleted.connect(self.file_deleted)
        self.bait2_sel_list_2.deleted.connect(self.file_deleted)
        self.bait2_nonsel_list_2.deleted.connect(self.file_deleted)

        # QThread for processes
        self.process = None
        self.set_interaction_state(False)
        self.verify_installation_btn.setEnabled(True)
        self.data = {}
        self.directory = '~/'
        self.jags_path = None
        self.r_path = None
        self.started = 0
        self.r = None
        self.fileio = f.fileio()
        self.stats = {}

    def initialize_folders(self, directory):
        self.fileio.create_new_folder(directory, "stat_maker_output")

    def pair_check(self, list1, list2):
        if list1.count() > 0 and list2.count() > 0:
            return True
        else:
            return False

    def either_check(self, list1, list2):
        if (list1.count() > 0 and list2.count() > 0) or (list1.count() == 0 and list2.count() == 0):
            return True
        else:
            return False

    def monitor_files(self):
        if self.pair_check(self.vec_sel_list, self.vec_nonsel_list) and \
                self.pair_check(self.bait1_sel_list, self.bait1_nonsel_list) and \
                self.pair_check(self.vec_sel_list_2, self.vec_nonsel_list_2):
            if self.either_check(self.bait2_sel_list, self.bait2_nonsel_list) and \
                    self.either_check(self.vec_sel_list_2, self.vec_nonsel_list_2) and \
                    self.either_check(self.bait1_sel_list_2, self.bait1_nonsel_list_2) and \
                    self.either_check(self.bait2_sel_list_2, self.bait2_nonsel_list_2):
                self.run_btn.setEnabled(True)
                self.overdisp_btn.setEnabled(True)
            else:
                self.run_btn.setEnabled(False)
                self.overdisp_btn.setEnabled(False)

    @QtCore.pyqtSlot()
    def on_overdisp_btn_clicked(self):
        QtWidgets.QMessageBox.information(self, "Statistics Engine Pending",
            "Overdispersion calculation (previously powered by R/JAGS) has been disabled in this build.\n\n"
            "It will be replaced with a DESeq2-based statistics pipeline in a future update.")

    def calculate_overdispersion(self):
        # R/JAGS-based overdispersion calculation is disabled pending replacement
        # with a DESeq2-based pipeline. See on_overdisp_btn_clicked.
        pass

    @QtCore.pyqtSlot()
    def on_verify_installation_btn_clicked(self):
        self.set_interaction_state(True)
        self.verify_installation_btn.setEnabled(False)
        self.verify_installation_btn.setText("Statistics Engine Not Yet Available")
        self.statusbar.showMessage("File selection is ready. Statistics computation (R/JAGS) is disabled pending "
                                   "replacement with a new pipeline.")

    def which(self, program):
        try:
            cmd = "which " + program
            ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = ps.communicate()[0]
            if len(output) > 0:
                return True
        except OSError:
            return False
        return False

    def get_path(self, program):
        try:
            cmd = "which " + program
            ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = ps.communicate()[0]
            if len(output) > 0:
                return output.rstrip().split()[0]
        except OSError:
            return None

    def file_deleted(self, path):
        for key in list(self.data.keys()):
            if self.data[key] == path:
                self.data.__delitem__(key)

    def fileDropped(self, list, url, key):
        if os.path.exists(url) and self.check_uniqueness(url):
            iconProvider = QtWidgets.QFileIconProvider()
            fileInfo = QtCore.QFileInfo(url)
            icon = iconProvider.icon(fileInfo)
            item = MyListWidgetItem()
            item.data = url
            item.setIcon(icon)
            item.setText(os.path.basename(os.path.normpath(url)))
            list.clear()
            list.addItem(item)
            self.data[key] = url
        self.monitor_files()

    def check_uniqueness(self, path):
        for url in self.data.values():
            if path == url:
                return False
        return True

    @QtCore.pyqtSlot()
    def on_quit_btn_clicked(self):
        app.quit()
        sys.exit()

    @QtCore.pyqtSlot()
    def on_folder_choice_btn_clicked(self):
        directory = str(QtWidgets.QFileDialog.getExistingDirectory(QtWidgets.QFileDialog(), "Locate Work Folder",
                                                               os.path.expanduser("~"),
                                                               QtWidgets.QFileDialog.ShowDirsOnly))
        self.directory = directory
        self.load_gene_summary_files()
        sys.excepthook = self.excepthook

    def excepthook(self, excType, excValue, tracebackobj):
        """
        Global function to catch unhandled exceptions.

        @param excType exception type
        @param excValue exception value
        @param tracebackobj traceback object
        """
        separator = '-' * 80
        logFile = self.directory + "/stat_maker_error.log"

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



    def load_gene_summary_files(self):
        self.initialize_folders(self.directory)
        try:
            dirlist = os.listdir(os.path.join(self.directory, 'gene_count_summary'))
            self.file_list.clear()
            for file in dirlist:
                if not re.match(r'^\.', file) and re.match(r'.+summary\.csv', file):
                    path = os.path.join(self.directory, 'gene_count_summary', file)
                    fileInfo = QtCore.QFileInfo(path)
                    iconProvider = QtWidgets.QFileIconProvider()
                    icon = iconProvider.icon(fileInfo)
                    item = MyListWidgetItem()
                    item.data = path
                    item.setIcon(icon)
                    item.setText(file)
                    self.file_list.addItem(item)
        except OSError:
            pass

    def vec_sel_fileDropped(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.vec_sel_list, path, 'Vector_Selected_1')

    def vec_nonsel_fileDropped(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.vec_nonsel_list, path, "Vector_Non-Selected_1")

    def bait1_sel_fileDropped(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.bait1_sel_list, path, "Bait1_Selected_1")

    def bait1_nonsel_fileDropped(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.bait1_nonsel_list, path, "Bait1_Non-Selected_1")

    def bait2_sel_fileDropped(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.bait2_sel_list, path, "Bait2_Selected_1")

    def bait2_nonsel_fileDropped(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.bait2_nonsel_list, path, "Bait2_Non-Selected_1")

    def vec_sel_fileDropped_2(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.vec_sel_list_2, path, "Vector_Selected_2")

    def vec_nonsel_fileDropped_2(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.vec_nonsel_list_2, path, "Vector_Non-Selected_2")

    def bait1_sel_fileDropped_2(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.bait1_sel_list_2, path, "Bait1_Selected_2")

    def bait1_nonsel_fileDropped_2(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.bait1_nonsel_list_2, path, "Bait1_Non-Selected_2")

    def bait2_sel_fileDropped_2(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.bait2_sel_list_2, path, "Bait2_Selected_2")

    def bait2_nonsel_fileDropped_2(self, path):
        path = str(path)
        if os.path.exists(path):
            self.fileDropped(self.bait2_nonsel_list_2, path, "Bait2_Non-Selected_2")

    def write_r_input(self):
        output = open(os.path.join(self.directory, 'r_input.params'), 'w')
        for key in sorted(self.data.keys()):
            output.write("%-25s = %s\n" % (key, self.data[key]))
        output.write("%-25s = %d\n" % ("Threshold", self.threshold_sbx.value()))
        output.close()

    def set_interaction_state(self, state):
        self.run_btn.setEnabled(state)
        self.overdisp_btn.setEnabled(state)
        self.file_list.setEnabled(state)
        self.vec_sel_list.setEnabled(state)
        self.vec_nonsel_list.setEnabled(state)
        self.bait1_sel_list.setEnabled(state)
        self.bait1_nonsel_list.setEnabled(state)
        self.bait2_sel_list.setEnabled(state)
        self.bait2_nonsel_list.setEnabled(state)

        self.vec_sel_list_2.setEnabled(state)
        self.vec_nonsel_list_2.setEnabled(state)
        self.bait1_sel_list_2.setEnabled(state)
        self.bait1_nonsel_list_2.setEnabled(state)
        self.bait2_sel_list_2.setEnabled(state)
        self.bait2_nonsel_list_2.setEnabled(state)
        self.folder_choice_btn.setEnabled(state)
        self.quit_btn.setEnabled(state)
        self.threshold_sbx.setEnabled(state)

    def _merge_dicts(*dict_args):
        """
        Given any number of dicts, shallow copy and merge into a new dict,
        precedence goes to key value pairs in latter dicts.
        """
        result = {}
        for dictionary in dict_args:
            result.update(dictionary)
        return result

    def get_junction_stats(self, key, name):
        gene_stat = {'frame_orf' : 0,
                     'upstream'  : 0,
                     'in_orf'    : 0,
                     'downstream': 0,
                     'in_frame'  : 0,
                     'backwards' : 0,
                     'intron'    : 0,
                     'total'     : 0,
                     'junctions' : {}
                     }
        try:
            gene = self.stats[key][name]
            for nm in gene.keys():
                if nm != 'stats':
                    gene_stat['junctions'][nm] = []
                    frame = 0
                    orf = 0
                    for j in gene[nm]:
                        if j.orf == 'in_orf' and j.frame == 'in_frame':
                            gene_stat['frame_orf'] += j.ppm
                        # Check orf
                        if j.orf == 'in_orf':
                            gene_stat['in_orf'] += j.ppm
                            orf = 2
                        elif j.orf == 'upstream':
                            gene_stat['upstream'] += j.ppm
                            orf = 1
                        elif j.orf == 'downstream':
                            gene_stat['downstream'] += j.ppm
                            orf = 3
                        # Check frame
                        if j.frame == 'in_frame':
                            gene_stat['in_frame'] += j.ppm
                            frame = 1
                        elif j.frame == 'backwards':
                            gene_stat['backwards'] += j.ppm
                            frame = 2
                        elif j.frame == 'intron':
                            gene_stat['intron'] += j.ppm
                            frame = 3
                        gene_stat['total'] += j.ppm
                        gene_stat['junctions'][nm].append([j.position,
                                                           j.query_start,
                                                           j.ppm, frame, orf])
        except KeyError:
            pass
        try:
            gene_stat['frame_orf'] = format(gene_stat['frame_orf'] * 100.0 / gene_stat['total'], ".1f")
            gene_stat['upstream'] = format(gene_stat['upstream'] * 100.0 / gene_stat['total'], ".1f")
            gene_stat['in_orf'] = format(gene_stat['in_orf'] * 100.0 / gene_stat['total'], ".1f")
            gene_stat['downstream'] = format(gene_stat['downstream'] * 100.0 / gene_stat['total'], ".1f")
            gene_stat['in_frame'] = format(gene_stat['in_frame'] * 100.0 / gene_stat['total'], ".1f")
            gene_stat['backwards'] = format(gene_stat['backwards'] * 100.0 / gene_stat['total'], ".1f")
            gene_stat['intron'] = format(gene_stat['intron'] * 100.0 / gene_stat['total'], ".1f")
        except ZeroDivisionError:
            gene_stat['frame_orf'] = '0'
            gene_stat['upstream'] = '0'
            gene_stat['in_orf'] = '0'
            gene_stat['downstream'] = '0'
            gene_stat['in_frame'] = '0'
            gene_stat['backwards'] = '0'
            gene_stat['intron'] = '0'
        return gene_stat

    def get_vector_gene_stats(self, key):
        filehandle = open(self.data[key])
        vec_genes = {}
        read = False
        for line in filehandle.readlines():
            if read:
                line_split = line.split(",")
                gene_name = line_split[1].lstrip().rstrip()
                vec_genes[gene_name] = self.get_junction_stats(key, gene_name)
                vec_genes[gene_name]['ppm'] = line_split[2].lstrip().rstrip()
            if re.match("^Chromosome", line):
                read = True
        filehandle.close()
        return vec_genes

    def runr(self):
        # R/JAGS-based statistics computation is disabled pending replacement
        # with a DESeq2-based pipeline. See on_run_btn_clicked.
        pass

    @QtCore.pyqtSlot()
    def on_run_btn_clicked(self):
        QtWidgets.QMessageBox.information(self, "Statistics Engine Pending",
            "Statistics computation (previously powered by R/JAGS) has been disabled in this build.\n\n"
            "It will be replaced with a DESeq2-based statistics pipeline in a future update.\n\n"
            "File selection and all other DEEPN tools are fully functional.")

def appExit():
    app.quit()
    sys.exit()

if __name__ == '__main__':
    form = Stat_Maker_Gui()
    form.show()
    app.aboutToQuit.connect(appExit)
    app.exec_()
