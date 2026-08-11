#!/usr/bin/env python3
import os
import sys

if '.app/Contents/MacOS' in sys.executable:
    os.chdir(os.path.join(os.path.dirname(sys.executable), '..', 'Resources'))

import re
import time
import traceback
import io
from time import localtime, strftime
import threading
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from collections import OrderedDict
import pickle
import xlsxwriter as xls

import functions.fileio_gui as f
import functions.plot as plot
import functions.structures as strt
import subprocess

app = QtWidgets.QApplication(sys.argv)
ui_path = os.path.join('ui', 'Query_Blast.ui')
if sys.platform == 'win32':
    ui_path = os.path.join('ui', 'Windows', 'Query_Blast.ui')
form_class, base_class = uic.loadUiType(ui_path)

main_directory = sys.argv[1]
list_name = sys.argv[2]

class QCustomTableWidgetItemInt(QtWidgets.QTableWidgetItem):
    def __init__(self, value):
        super(QCustomTableWidgetItemInt, self).__init__('%d' % value)

    def __lt__(self, other):
        if isinstance(other, QCustomTableWidgetItemInt):
            self_data_value = int(str(self.data(QtCore.Qt.EditRole)))
            other_data_value = int(str(other.data(QtCore.Qt.EditRole)))
            return self_data_value < other_data_value
        else:
            return QtWidgets.QTableWidgetItem.__lt__(self, other)

class QCustomTableWidgetItemFloat(QtWidgets.QTableWidgetItem):
    def __init__(self, value):
        super(QCustomTableWidgetItemFloat, self).__init__('%.3f' % value)

    def __lt__(self, other):
        if isinstance(other, QCustomTableWidgetItemFloat):
            self_data_value = float(str(self.data(QtCore.Qt.EditRole)))
            other_data_value = float(str(other.data(QtCore.Qt.EditRole)))
            return self_data_value < other_data_value
        else:
            return QtWidgets.QTableWidgetItem.__lt__(self, other)

class Query_Blast_Gui(QtWidgets.QMainWindow, form_class):
    def __init__(self, folder, directory, list_file, *args):
        super(Query_Blast_Gui, self).__init__(*args)
        self.setupUi(self)
        self.clipboard = ""
        self.input_folder = folder
        self.directory = directory
        self.list_file = list_file
        self.updating = False
        self.fileio = f.fileio()
        # Selections from GUI
        self.selected_gene_name = ''
        self.selected_accession_number = ''
        self.selected_dataset_names = ['-- Select --', '-- Select --', '-- Select --']
        # Generic Table Names
        self.table = None
        self.table_filtered = None
        # Cache for storing data
        self.gene_info_list = OrderedDict()
        self.accession_number_list = OrderedDict()
        self.datasets_cache = OrderedDict()
        self.sequences_info_list = {}
        # Methods
        self._initialize_ui_elements()
        # Results Store
        self.results = [[], [], []]
        self.filtered_results = [[], [], []]

        self.graph1 = plot.QBPlot('Base Pair', 'Junction Count (ppm)')
        self.tabWidget_1.addTab(self.graph1, "Plot")
        self.graph1.initialize_plot()

        self.graph2 = plot.QBPlot('Base Pair', 'Junction Count (ppm)')
        self.tabWidget_2.addTab(self.graph2, "Plot")
        self.graph2.initialize_plot()

        self.graph3 = plot.QBPlot('Base Pair', 'Junction Count (ppm)')
        self.tabWidget_3.addTab(self.graph3, "Plot")
        self.graph3.initialize_plot()

    def monitor_clipboard(self):
        while 1:
            self.clipboard = str(subprocess.check_output('pbpaste'))
            self.clipboard = self.clipboard.rstrip().lstrip()
            if self.clipboard in self.gene_info_list.keys():
                if str(self.search_bar.text()).rstrip().lstrip() != self.clipboard:
                    print("Updating Results for Gene... ", self.clipboard)
                    sys.stdout.flush()
                    self.search_bar.setText(self.clipboard)
                    self.old_clipboard = self.clipboard
                    self.search_bar.editingFinished.emit()
                else:
                    subprocess.check_output('pbcopy </dev/null', shell=True)
            time.sleep(0.5)

    def _initialize_ui_elements(self):
        threading.Thread(target=self.monitor_clipboard, daemon=True).start()
        self.get_sequences()
        self.populate_gene_suggestions()
        self.blast_dataset_ddl_1.addItem('-- Select --')
        self.blast_dataset_ddl_2.addItem('-- Select --')
        self.blast_dataset_ddl_3.addItem('-- Select --')
        for fi in self.fileio.get_file_list(self.directory, self.input_folder, '.bqp'):
            self.blast_dataset_ddl_1.addItem(fi)
            self.blast_dataset_ddl_2.addItem(fi)
            self.blast_dataset_ddl_3.addItem(fi)
        self.blast_dataset_ddl_1.currentIndexChanged.connect(self._dataset_selected)
        self.blast_dataset_ddl_2.currentIndexChanged.connect(self._dataset_selected)
        self.blast_dataset_ddl_3.currentIndexChanged.connect(self._dataset_selected)
        self.disable_ui_elements()

    def populate_gene_suggestions(self):
        self.completer = QtWidgets.QCompleter(list(self.accession_number_list.keys()) + list(self.gene_info_list.keys()), self)
        self.completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        self.search_bar.setCompleter(self.completer)
        self.search_bar.editingFinished.connect(self.search_changed)
        self.accession_list.currentIndexChanged.connect(self.accession_changed)

    def search_changed(self):
        self.accession_list.clear()
        self.search_text = str(self.search_bar.text()).rstrip().lstrip()
        if self.search_text in self.gene_info_list.keys():
            self.selected_gene_name = self.search_text
            items = [self.accession_list.itemText(i) for i in range(self.accession_list.count())]
            for accession in self.gene_info_list[self.selected_gene_name]:
                if accession not in items:
                    self.accession_list.addItem(accession)
            self.accession_list.setCurrentIndex(0)
            self.enable_ui_elements()
        elif self.search_text in self.accession_number_list.keys():
            self.selected_gene_name = self.accession_number_list[self.search_text]
            for accession in self.gene_info_list[self.selected_gene_name]:
                self.accession_list.addItem(accession)
            self.search_bar.setText(self.selected_gene_name)
            self.accession_list.setCurrentIndex(0)
            self.enable_ui_elements()
        else:
            self.disable_ui_elements()
            self.status_bar.showMessage("Gene Not Found!", 5000)
            self.sequence_browser.clear()
            self.graph1.clear_plot()
            self.graph2.clear_plot()
            self.graph3.clear_plot()
            if self.table:
                self.table.clearContents()
                self.table.setRowCount(0)

            if self.table_filtered:
                self.table_filtered.clearContents()
                self.table_filtered.setRowCount(0)

        if self.selected_gene_name != '' and str(self.accession_list.currentText()) != '':
            self._dataset_selected()
            _seq = self.sequence_format_html(self.sequences_info_list[str(self.accession_list.currentText())])
            self.sequence_browser.setHtml(_seq)

    def disable_ui_elements(self):
        self.sequence_browser.setEnabled(False)
        self.blast_dataset_ddl_1.setEnabled(False)
        self.blast_dataset_ddl_2.setEnabled(False)
        self.blast_dataset_ddl_3.setEnabled(False)
        self.tabWidget_1.setEnabled(False)
        self.tabWidget_2.setEnabled(False)
        self.tabWidget_3.setEnabled(False)

    def enable_ui_elements(self):
        self.sequence_browser.setEnabled(True)
        self.blast_dataset_ddl_1.setEnabled(True)
        self.blast_dataset_ddl_2.setEnabled(True)
        self.blast_dataset_ddl_3.setEnabled(True)
        self.tabWidget_1.setEnabled(True)
        self.tabWidget_2.setEnabled(True)
        self.tabWidget_3.setEnabled(True)

    def accession_changed(self):
        if self.selected_gene_name != '' and str(self.accession_list.currentText()) != '':
            self.selected_accession_number = str(self.accession_list.currentText())
            if self.selected_accession_number != '':
                _seq = self.sequence_format_html(self.sequences_info_list[self.selected_accession_number])
                self.sequence_browser.setHtml(_seq)
            self._dataset_selected()


    def sequence_format_html(self, selected_gene):
        sequence = selected_gene['mRNA']
        start = selected_gene['orf_start'] - 1
        stop = selected_gene['orf_stop']
        html_string = '<html><body>> %s (%s) on chromosome %s (%s) ORF Start: %d ORF Stop: %d<br><font ' \
                      'color="silver">' % (selected_gene['nm_number'], selected_gene['gene_name'],
                                           selected_gene['chromosome'], selected_gene['intron'],
                                           selected_gene['orf_start'], selected_gene['orf_stop'])
        html_string += sequence[:start] + '</font>'
        html_string += sequence[start:stop]
        html_string += '<font color="silver">' + sequence[stop:] + '</font>'
        html_string += '</body></html>'
        return html_string

    def get_sequences(self):
        self.gene_info_list = OrderedDict()
        self.accession_number_list = OrderedDict()
        fh = open(os.path.join('lists', self.list_file))
        self.sequences_info_list = {}
        for line in fh.readlines():
            split = line.split()
            self.accession_number_list[split[0]] = split[1]
            if split[1] not in self.gene_info_list:
                self.gene_info_list[split[1]] = []
            self.gene_info_list[split[1]].append("%s" % split[0])
            self.sequences_info_list[split[0]] = {'gene_name' : split[1],
                                                  'nm_number' : split[0],
                                                  'orf_start' : int(split[6]) + 1,
                                                  'orf_stop'  : int(split[7]),
                                                  'mRNA'      : split[9],
                                                  'intron'    : split[8],
                                                  'chromosome': split[2]
                                                  }

    def _initialize_results_store(self):
        self.results = [[], [], []]
        self.filtered_results = [[], [], []]
        self.results[0].append(('Position', "# Junctions", "Query Start", "CDS Status", "Reading Frame"))
        self.filtered_results[0].append(('Position', "# Junctions"))
        self.results[1].append(('Position', "# Junctions", "Query Start", "CDS Status", "Reading Frame"))
        self.filtered_results[1].append(('Position', "# Junctions"))
        self.results[2].append(('Position', "# Junctions", "Query Start", "CDS Status", "Reading Frame"))
        self.filtered_results[2].append(('Position', "# Junctions"))

    def _get_sender_index_name(self, sender):
        index = 0
        sender_objectname = sender.objectName()
        try:
            index = int(sender_objectname[-1:]) - 1
        except ValueError:
            pass
        dataset_name = '-- Select --'
        try:
            dataset_name = str(self.sender().currentText())
        except AttributeError:
            pass
        return index, dataset_name

    def make_condensed_list(self, inlist): #What does this function do???
        count = 0
        roster = set()
        return_list = []
        for item in inlist:
            if item[0] not in roster:
                roster.add(item[0])

        for item in roster:
            for thing in inlist:
                if thing[0] == item:
                    count += thing[2]
            return_list.append((item, count))
            count = 0
        return_list.sort()
        return return_list

    def get_filtered_data(self, junctions1, junctions2, table, index, filtered_results):
        freq1 = []
        freq2 = []
        for j1 in junctions1:
            freq1.append((j1.position, j1.query_start, j1.ppm))
        freq1B = self.make_condensed_list(freq1)
        for j2 in junctions2:
            freq2.append((j2.position, j2.query_start, j2.ppm))
        freq2B = self.make_condensed_list(freq2)
        for item2 in freq2B:
            if item2[0] in [item[0] for item in freq1B]:
                filtered_results[index].append((item2[0], round(item2[1], 3)))
        table.clearContents()
        table.setSortingEnabled(False)
        table.setRowCount(len(filtered_results[index]) - 1)
        for i, row in enumerate(filtered_results[index]):
            if i != 0:
                table.setItem(i-1, 0, QCustomTableWidgetItemInt(row[0]))
                table.setItem(i-1, 1, QCustomTableWidgetItemFloat(row[1]))
        table.setSortingEnabled(True)
        table.sortItems(1, 1)

    def load_dataset(self, dataset_name):
        self.datasets_cache[dataset_name] = pickle.load(open(os.path.join(self.directory, self.input_folder,
                                                                           dataset_name), 'rb'))

        self.enable_ui_elements()
        self._dataset_selected()

    def _dataset_selected(self):
        index = 0
        dataset_name = self.selected_dataset_names[index]
        self._initialize_results_store()

        if isinstance(self.sender(), QtWidgets.QComboBox):
            if self.sender().objectName() != 'accession_list':
                index, dataset_name = self._get_sender_index_name(self.sender())
                self.selected_dataset_names[index] = dataset_name
        else:
            pass

        if dataset_name != '-- Select --':
            try:
                self.datasets_cache[dataset_name]
            except KeyError:
                self.disable_ui_elements()
                print(dataset_name)
                threading.Thread(target=self.load_dataset, args=(dataset_name,), daemon=True).start()

            count = 1
            for name in self.selected_dataset_names:
                if name != '-- Select --':
                    dataset = self.datasets_cache[name]
                    try:
                        exec("%s = %s" % ("self.table", "self.results_lst_" + str(count)))
                        exec("%s = %s" % ("self.table_filtered", "self.results_filtered_lst_" + str(count)))
                        exec("%s = %s" % ("self.graph", "self.graph" + str(count)))
                    except AttributeError:
                        pass

                    self.table.clearContents()
                    junction_data = dataset[self.selected_gene_name][self.selected_accession_number]
                    self.table.setSortingEnabled(False)
                    self.table.setRowCount(len(junction_data))
                    for row, values in enumerate(junction_data):
                        j = junction_data[row]
                        self.table.setItem(row, 0, QCustomTableWidgetItemInt(j.position))
                        self.table.setItem(row, 1, QCustomTableWidgetItemFloat(round(j.ppm, 3)))
                        self.table.setItem(row, 2, QCustomTableWidgetItemInt(j.query_start))
                        self.table.setItem(row, 3, QtWidgets.QTableWidgetItem(j.orf.replace("_", " ").title()))
                        self.table.setItem(row, 4, QtWidgets.QTableWidgetItem(j.frame.replace("_", " ").title()))
                        data = (str(j.position), str(round(j.ppm, 3)), str(j.query_start), j.orf, j.frame)
                        self.results[count - 1].append(data)
                    self.table.setSortingEnabled(True)
                    self.table.sortItems(1, 1)
                    self.graph.plot(self.results[count - 1],
                                    self.sequences_info_list[str(self.accession_list.currentText())]['orf_start'],
                                    self.sequences_info_list[str(self.accession_list.currentText())]['orf_stop'])
                    try:
                        junction_data1 = self.datasets_cache[self.selected_dataset_names[0]][self.selected_gene_name][self.selected_accession_number]
                        junction_data2 = dataset[self.selected_gene_name][self.selected_accession_number]
                        self.get_filtered_data(junction_data1, junction_data2,
                                               self.table_filtered, count - 1,
                                               self.filtered_results)
                    except KeyError:
                        pass

                    count += 1

    @QtCore.pyqtSlot()
    def on_save_btn_clicked(self):
        self.save_query_results(self.results, self.filtered_results, self.selected_gene_name,
                                self.selected_dataset_names, self.directory)
        self.status_bar.showMessage("Saved to file %s" % os.path.join('query_results', '%s_%s.xlsx' %
                                    (self.selected_gene_name, strftime("%Y-%m-%d_%Hh%Mm%Ss",
                                                                       localtime()))),
                                    5000)

    def save_query_results(self, results, filtered_results, file_name, dataset_names, directory, path='query_results'):
        self.fileio.create_new_folder(directory, path)
        workbook = xls.Workbook(os.path.join(directory, path, '%s_%s.xlsx' %
                                             (file_name, strftime("%Y-%m-%d_%Hh%Mm%Ss", localtime()))))
        for index, result in enumerate(results):
            worksheet = workbook.add_worksheet(dataset_names[index][:20] + '..._' + str(index))
            for row, row_item in enumerate(result):
                for col, item in enumerate(row_item):
                    worksheet.write(row, col, item)

            for row, row_item in enumerate(filtered_results[index]):
                for col, item in enumerate(row_item):
                    worksheet.write(row, col+10, item)
        workbook.close()


def appExit():
    app.quit()
    sys.exit()

if __name__ == '__main__':
    input_folder = 'blast_results_query'
    form = Query_Blast_Gui(input_folder, main_directory, list_name)
    form.show()
    app.aboutToQuit.connect(appExit)
    app.exec_()
