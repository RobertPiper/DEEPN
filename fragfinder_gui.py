#!/usr/bin/env python3
import os
import sys

if '.app/Contents/MacOS' in sys.executable:
    os.chdir(os.path.join(os.path.dirname(sys.executable), '..', 'Resources'))

import pickle
import time
from collections import OrderedDict
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui, QtWidgets, uic

import functions.fileio_gui as f
import functions.plot as plot
import functions.process as process

app = QtWidgets.QApplication(sys.argv)
ui_path = os.path.join('ui', 'FragFinder.ui')
if sys.platform == 'win32':
    ui_path = os.path.join('ui', 'Windows', 'FragFinder.ui')
form_class, base_class = uic.loadUiType(ui_path)
pg.setConfigOption('background', (236, 236, 236))
pg.setConfigOption('foreground', 'k')

main_directory = sys.argv[1]
gene_list_file = sys.argv[2]
combined = int(sys.argv[3])

RD_LEGEND_ENTRIES = [('#1A64B0', 'In ORF'), ('#777777', 'Not in CDS'), ('#D13A26', 'Start/Stop of CDS')]
JUNCTION_LEGEND_ENTRIES = [('#1A64B0', 'In ORF'), ('#1CC5FF', 'Upstream'),
                           ('#777777', 'Downstream / Out of Frame'), ('#D13A26', 'Start/Stop of CDS')]


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


class GetReadList_Thread(QtCore.QThread):
    GeneListFinished = QtCore.pyqtSignal()

    def __init__(self, directory, input_folder, dataset_name, selected_gene, parent):
        QtCore.QThread.__init__(self, parent)
        self.directory = directory
        self.input_folder = input_folder
        self.dataset_name = dataset_name
        self.selected_gene = selected_gene
        self.parent = parent

    def run(self):
        filehandle = open(os.path.join(self.directory, 'gene_count_indices', self.dataset_name[:-4],
                                       self.selected_gene['chromosome'] + '.bin'), 'r')
        self.parent.list_of_reads = []
        for line in filehandle.readlines():
            split = line.strip().split(':')
            if int(split[0]) >= int(self.selected_gene['chromosome_start']) and int(split[0]) <= int(self.selected_gene['chromosome_stop']):
                self.parent.list_of_reads.append(split[1])
        self.GeneListFinished.emit()


class CalculateDepth_Thread(QtCore.QThread):
    CalculateDepthFinished = QtCore.pyqtSignal()

    def __init__(self, parent):
        QtCore.QThread.__init__(self, parent)
        self.process = process.process()
        self.parent = parent

    def make_list_of_mRNA(self):
        # .prn mRNA sequences are stored lowercase, but real reads (and
        # reverse_complement()'s output, which uppercases internally) are
        # uppercase - without this .upper(), the forward-strand comparison in
        # run() below silently never matches, undercounting read depth to
        # whatever the reverse-complement side alone picks up.
        sequence = []
        for position in range(0, len(self.parent.selected_gene['mRNA']), self.parent.interval):
            seq = self.parent.selected_gene['mRNA'][position:position + 20].upper()
            sequence.append((position, seq))
        return sequence

    def run(self):
        sequences = self.make_list_of_mRNA()
        self.parent.results = []
        for (position, item) in sequences:
            item_reverse = self.process.reverse_complement(item)
            count = 0
            for read in self.parent.list_of_reads:
                if (str(item) in str(read)) or (str(item_reverse) in str(read)):
                    count += 1
            orfness = "Not in CDS"
            if int(position) >= int(self.parent.selected_gene['orf_start']) and int(position) <= int(self.parent.selected_gene['orf_stop']):
                orfness = "In ORF"
            self.parent.results.append((position + 1, count, orfness))
        self.CalculateDepthFinished.emit()


class ExtrapolateJunctionThread(QtCore.QThread):
    """Re-scans a dataset's raw blast.txt / junctions.txt (rather than the
    aggregated .bqp) to find every individual read supporting one exact
    junction (position + query_start), reapplying the same hit-quality
    filter Junction Make used, and returns the longest one.

    Progress is reported by both percent and elapsed time, not percent
    alone: on a slow or throttled volume, minutes can pass between whole
    percentage points, and a dialog that only moves on percent change
    looks frozen even while it's actively scanning. Forcing an update at
    least once a second (with the actual MB scanned) keeps it visibly
    live regardless of how slow the underlying I/O is."""
    progress = QtCore.pyqtSignal(int, str)
    finished_result = QtCore.pyqtSignal(object, object, object, object, str)

    def __init__(self, blast_txt_path, junctions_txt_path, direction, accession, position, query_start, parent=None):
        QtCore.QThread.__init__(self, parent)
        self.blast_txt_path = blast_txt_path
        self.junctions_txt_path = junctions_txt_path
        self.direction = direction
        self.accession = accession
        self.position = position
        self.query_start = query_start
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def _emit_progress(self, pct, bytes_read, total_bytes, label):
        self.progress.emit(pct, "%s... %.0f MB / %.0f MB" % (label, bytes_read / 1e6, total_bytes / 1e6))

    def run(self):
        matching_read_ids = set()
        backwards_by_read = {}
        previous_bitscore = 0
        collect_results = 'n'

        file_size = os.path.getsize(self.blast_txt_path) or 1
        last_pct = -1
        last_emit_time = time.time()
        blast_fh = open(self.blast_txt_path, 'r')
        while True:
            line = blast_fh.readline()
            if not line:
                break
            if self._cancelled:
                blast_fh.close()
                self.finished_result.emit(None, None, None, None, self.direction)
                return
            split = line.split()
            if not split:
                continue
            if "BLASTN" in line:
                previous_bitscore = 0
            elif "hits" in line and int(split[1]) < 100:
                collect_results = 'y'
            elif split[0] != '#' and collect_results == 'y' and float(split[2]) > 98 and \
                            float(split[11]) > 50.0 and float(split[11]) > previous_bitscore:
                previous_bitscore = float(split[11]) * 0.98
                if split[1] == self.accession and int(split[8]) == self.position and int(split[6]) == self.query_start:
                    matching_read_ids.add(split[0])
                    backwards_by_read[split[0]] = (int(split[9]) - int(split[8])) < 0
            else:
                collect_results = 'n'

            pos = blast_fh.tell()
            pct = min(49, int(pos * 49.0 / file_size))
            now = time.time()
            if pct != last_pct or (now - last_emit_time) >= 1.0:
                self._emit_progress(pct, pos, file_size, "Scanning BLAST results")
                last_pct = pct
                last_emit_time = now
        blast_fh.close()

        if not matching_read_ids:
            self.finished_result.emit(None, None, None, None, self.direction)
            return

        best_read_id = None
        best_seq = ''
        best_full_read = None
        file_size2 = os.path.getsize(self.junctions_txt_path) or 1
        last_pct = 49
        last_emit_time = time.time()
        junctions_fh = open(self.junctions_txt_path, 'r')
        while True:
            line = junctions_fh.readline()
            if not line:
                break
            if self._cancelled:
                junctions_fh.close()
                self.finished_result.emit(None, None, None, None, self.direction)
                return
            split = line.split()
            if len(split) > 5 and split[0] in matching_read_ids:
                seq = split[5]
                if len(seq) > len(best_seq):
                    best_seq = seq
                    best_full_read = split[4]
                    best_read_id = split[0]
            pos = junctions_fh.tell()
            pct = 50 + min(49, int(pos * 49.0 / file_size2))
            now = time.time()
            if pct != last_pct or (now - last_emit_time) >= 1.0:
                self._emit_progress(pct, pos, file_size2, "Scanning junction reads")
                last_pct = pct
                last_emit_time = now
        junctions_fh.close()

        self.progress.emit(100, "Done.")
        is_backwards = backwards_by_read.get(best_read_id, False)
        self.finished_result.emit(best_read_id, best_seq, best_full_read, is_backwards, self.direction)


class FragFinder_Gui(QtWidgets.QMainWindow, form_class):
    def __init__(self, directory, gene_list_file, combined, *args):
        super(FragFinder_Gui, self).__init__(*args)
        self.setupUi(self)
        self.directory = directory
        self.gene_list_file = gene_list_file
        self.combined = combined
        self.fileio = f.fileio()

        self.selected_gene_name = ''
        self.selected_accession_number = ''
        self.selected_gene = {}
        self.gene_list = OrderedDict()
        self.gene_info_list = OrderedDict()
        self.accession_number_list = OrderedDict()

        self.rd_base = None
        self.fivep_base = None
        self.threep_base = None
        self.fivep_data = None
        self.threep_data = None
        # The actual prefix found on the loaded filename ('5p_', '3p_', or ''
        # for data processed before v4's prefix convention existed) - tracked
        # separately per panel since a legacy unprefixed file can be loaded
        # into either panel, and Extrapolate Junction needs this exact prefix
        # to find the matching raw blast.txt/junctions.txt files.
        self.fivep_file_prefix = ''
        self.threep_file_prefix = ''
        self.table_row_data = []
        self.selected_junction = None

        self.interval = 100
        self.list_of_reads = []
        self.results = []

        self.rd_plot = plot.RDPlot("Base Pair", "# Hits", show_legend=False)
        self.fivep_plot = plot.JunctionPlot("Base Pair", "Junction Count (ppm)", show_legend=False)
        self.threep_plot = plot.JunctionPlot("Base Pair", "Junction Count (ppm)", show_legend=False)
        self.plot_layout_rd.addWidget(self.rd_plot)
        self.plot_layout_rd.addWidget(plot.make_legend_widget(RD_LEGEND_ENTRIES))
        self.plot_layout_5p.addWidget(self.fivep_plot)
        self.plot_layout_5p.addWidget(plot.make_legend_widget(JUNCTION_LEGEND_ENTRIES))
        self.plot_layout_3p.addWidget(self.threep_plot)
        self.plot_layout_3p.addWidget(plot.make_legend_widget(JUNCTION_LEGEND_ENTRIES))
        self.fivep_plot.setXLink(self.rd_plot)
        self.threep_plot.setXLink(self.rd_plot)
        # Vertical stretch factors alone don't stop these plots from also
        # absorbing extra window height, since they're Expanding widgets -
        # Qt's box layout treats stretch as a weight between competing
        # Expanding items, not a hard 0-or-nothing switch. Pin each plot to
        # its natural height (captured after the first real layout pass) and
        # switch it to a Fixed vertical policy, so any window growth goes
        # entirely to sequence_browser instead.
        QtCore.QTimer.singleShot(0, self._lock_plot_heights)

        self.gene_list, self.gene_info_list, self.accession_number_list = self.get_gene_and_accession_lists(self.gene_list_file)
        self.populate_accession_suggestions()

        # rd_load_btn/fivep_load_btn/threep_load_btn/extrapolate_btn/copy_btn are
        # deliberately NOT connected here - their handler names already match
        # PyQt's on_<objectName>_<signalName> auto-connect convention, which
        # setupUi() wires up on its own; connecting them again here fired each
        # handler twice per click (e.g. two Extrapolate Junction progress dialogs).
        self.junction_table.itemSelectionChanged.connect(self.on_junction_table_selection_changed)
        self.rd_interval_spin.setValue(self.interval)
        self.rd_interval_spin.valueChanged.connect(self.on_rd_interval_spin_changed)

        self.read_list_thread = GetReadList_Thread(self.directory, None, None, self.selected_gene, self)
        self.read_list_thread.GeneListFinished.connect(self.on_rd_read_list_finished)
        self.calculate_depth_thread = CalculateDepth_Thread(self)
        self.calculate_depth_thread.CalculateDepthFinished.connect(self.on_rd_depth_calculated)

    def _lock_plot_heights(self):
        for p in (self.rd_plot, self.fivep_plot, self.threep_plot):
            p.setFixedHeight(p.height())
            p.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)

    # ---- gene search (same pattern as Read Depth v2 / Blast Query) ----

    def get_gene_and_accession_lists(self, gene_list_file):
        fh = open(os.path.join('lists', gene_list_file), 'r')
        gene_list = OrderedDict()
        gene_info_list = OrderedDict()
        accession_number_list = OrderedDict()
        for line in fh.readlines():
            split = line.split()
            accession = split[0]
            gene_name = split[1]
            gene_list[accession] = {'nm_number': accession, 'gene_name': gene_name, 'orf_start': int(split[6]) + 1,
                                   'orf_stop': int(split[7]), 'mRNA': split[9], 'intron': split[8],
                                   'chromosome': split[2], 'chromosome_start': split[4], 'chromosome_stop': split[5]}
            accession_number_list[accession] = gene_name
            if gene_name not in gene_info_list:
                gene_info_list[gene_name] = []
            gene_info_list[gene_name].append(accession)
        return gene_list, gene_info_list, accession_number_list

    def populate_accession_suggestions(self):
        self.completer = QtWidgets.QCompleter(list(self.gene_info_list.keys()) + list(self.accession_number_list.keys()), self)
        self.completer.setCompletionMode(QtWidgets.QCompleter.PopupCompletion)
        self.completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        self.search_bar.setCompleter(self.completer)
        self.search_bar.editingFinished.connect(self.search_changed)
        self.accession_list.currentIndexChanged.connect(self.accession_changed)

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

    @QtCore.pyqtSlot()
    def search_changed(self):
        self.accession_list.clear()
        search_text = str(self.search_bar.text()).strip()
        if search_text in self.gene_info_list.keys():
            self.selected_gene_name = search_text
            for accession in self.gene_info_list[self.selected_gene_name]:
                self.accession_list.addItem(accession)
            self.accession_list.setCurrentIndex(0)
        elif search_text in self.accession_number_list.keys():
            self.selected_gene_name = self.accession_number_list[search_text]
            self.search_bar.setText(self.selected_gene_name)
            for accession in self.gene_info_list[self.selected_gene_name]:
                self.accession_list.addItem(accession)
            index = self.accession_list.findText(search_text)
            self.accession_list.setCurrentIndex(index if index >= 0 else 0)
        else:
            self.status_bar.showMessage("Gene Not Found!", 5000)
            self.selected_gene_name = ''
            self.selected_accession_number = ''
            self.selected_gene = {}
            self.sequence_browser.clear()
            self.update_plots_for_selected_gene()

    @QtCore.pyqtSlot()
    def accession_changed(self):
        self.selected_accession_number = str(self.accession_list.currentText())
        if self.selected_accession_number == '':
            return
        try:
            self.selected_gene = self.gene_list[self.selected_accession_number]
        except KeyError:
            self.selected_gene = {}
            return
        self.sequence_browser.setHtml(self.sequence_format_html(self.selected_gene))
        self.update_plots_for_selected_gene()

    # ---- file discovery / load-and-prefill ----

    def _sam_folder(self):
        return 'sam_files' if self.combined == 1 else 'mapped_sam_files'

    @QtCore.pyqtSlot()
    def on_rd_interval_spin_changed(self):
        self.interval = self.rd_interval_spin.value()
        if self.rd_base is not None and self.selected_gene and self.selected_accession_number:
            self.rd_load_btn.setEnabled(False)
            self.status_bar.showMessage("Calculating read depth...")
            self.calculate_depth_thread.start()

    @QtCore.pyqtSlot()
    def on_rd_load_btn_clicked(self):
        folder = os.path.join(self.directory, self._sam_folder())
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Locate Read Depth .sam file", folder, "SAM files (*.sam)")
        if not path:
            return
        filename = os.path.basename(path)
        base = filename[:-4] if filename.endswith('.sam') else filename
        self._load_rd(filename, base)
        self._prefill_companions(base, skip='rd')

    @QtCore.pyqtSlot()
    def on_fivep_load_btn_clicked(self):
        self._load_bqp_via_dialog('5p')

    @QtCore.pyqtSlot()
    def on_threep_load_btn_clicked(self):
        self._load_bqp_via_dialog('3p')

    @staticmethod
    def _split_bqp_prefix(filename):
        """Splits a .bqp filename into (prefix, base). Recognizes the v4
        5p_/3p_ convention; anything else (including data processed before
        that convention existed) is treated as unprefixed, with the whole
        filename as the base."""
        for candidate_prefix in ('5p_', '3p_'):
            if filename.startswith(candidate_prefix):
                return candidate_prefix, filename[len(candidate_prefix):-4]
        return '', filename[:-4]

    def _load_bqp_via_dialog(self, direction):
        label = "5'" if direction == '5p' else "3'"
        folder = os.path.join(self.directory, 'blast_results_query')
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Locate %s Junction .bqp file" % label, folder,
                                                         "All .bqp files (*.bqp);;All files (*)")
        if not path:
            return
        filename = os.path.basename(path)
        file_prefix, base = self._split_bqp_prefix(filename)
        self._load_bqp(direction, filename, base, file_prefix)
        self._prefill_companions(base, skip=direction)

    def _load_rd(self, filename, base):
        self.rd_base = base
        self.rd_filename_lbl.setText(': ' + filename)
        self.update_plots_for_selected_gene()

    def _load_bqp(self, direction, filename, base, file_prefix):
        data = pickle.load(open(os.path.join(self.directory, 'blast_results_query', filename), 'rb'))
        if direction == '5p':
            self.fivep_base = base
            self.fivep_data = data
            self.fivep_file_prefix = file_prefix
            self.fivep_filename_lbl.setText(': ' + filename)
        else:
            self.threep_base = base
            self.threep_data = data
            self.threep_file_prefix = file_prefix
            self.threep_filename_lbl.setText(': ' + filename)
        self.update_plots_for_selected_gene()

    def _prefill_companions(self, base, skip):
        if skip != 'rd':
            self._auto_or_prompt_rd(base)
        if skip != '5p':
            self._auto_or_prompt_bqp('5p', base)
        if skip != '3p':
            self._auto_or_prompt_bqp('3p', base)

    def _auto_or_prompt_rd(self, base):
        folder = os.path.join(self.directory, self._sam_folder())
        candidate = base + '.sam'
        if os.path.exists(os.path.join(folder, candidate)):
            self._load_rd(candidate, base)
            return
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Locate matching Read Depth .sam file (Cancel to skip)",
                                                         folder, "SAM files (*.sam);;All files (*)")
        if path:
            filename = os.path.basename(path)
            self._load_rd(filename, filename[:-4] if filename.endswith('.sam') else filename)

    def _auto_or_prompt_bqp(self, direction, base):
        label = "5'" if direction == '5p' else "3'"
        prefix = direction + '_'
        folder = os.path.join(self.directory, 'blast_results_query')
        candidate = prefix + base + '.bqp'
        if os.path.exists(os.path.join(folder, candidate)):
            self._load_bqp(direction, candidate, base, prefix)
            return
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Locate matching %s Junction .bqp file (Cancel to skip)" % label,
                                                         folder, "%s junction files (%s*.bqp);;All .bqp files (*.bqp);;All files (*)" %
                                                         (label, prefix))
        if path:
            filename = os.path.basename(path)
            file_prefix, picked_base = self._split_bqp_prefix(filename)
            self._load_bqp(direction, filename, picked_base, file_prefix)

    # ---- plots + junction table ----

    def update_plots_for_selected_gene(self):
        self.junction_table.setSortingEnabled(False)
        self.junction_table.clearContents()
        self.junction_table.setRowCount(0)
        self.junction_table.setSortingEnabled(True)
        self.table_row_data = []
        self.extrapolate_btn.setEnabled(False)
        self.selected_junction = None
        self.rd_plot.clear_marker()
        self.fivep_plot.clear_marker()
        self.threep_plot.clear_marker()

        if not self.selected_gene or not self.selected_accession_number:
            self.rd_plot.clear_plot()
            self.fivep_plot.clear_plot()
            self.threep_plot.clear_plot()
            return

        start = self.selected_gene['orf_start']
        stop = self.selected_gene['orf_stop']

        if self.rd_base is not None:
            self.rd_load_btn.setEnabled(False)
            self.status_bar.showMessage("Calculating read depth...")
            self.read_list_thread.dataset_name = self.rd_base + '.sam'
            self.read_list_thread.selected_gene = self.selected_gene
            self.read_list_thread.start()
        else:
            self.rd_plot.clear_plot()

        rows = []
        if self.fivep_data is not None:
            junctions = self.fivep_data.get(self.selected_gene_name, {}).get(self.selected_accession_number, [])
            self.fivep_plot.plot(junctions, start, stop)
            for j in junctions:
                rows.append(('5p', j))
        else:
            self.fivep_plot.clear_plot()

        if self.threep_data is not None:
            junctions = self.threep_data.get(self.selected_gene_name, {}).get(self.selected_accession_number, [])
            self.threep_plot.plot(junctions, start, stop)
            for j in junctions:
                rows.append(('3p', j))
        else:
            self.threep_plot.clear_plot()

        self._populate_junction_table(rows)

    def _populate_junction_table(self, rows):
        self.table_row_data = rows
        self.junction_table.setSortingEnabled(False)
        self.junction_table.setRowCount(len(rows))
        for i, (direction, j) in enumerate(rows):
            label = "5'" if direction == '5p' else "3'"
            position_item = QCustomTableWidgetItemInt(j.position)
            position_item.setData(QtCore.Qt.UserRole, i)
            self.junction_table.setItem(i, 0, position_item)
            self.junction_table.setItem(i, 1, QtWidgets.QTableWidgetItem(label))
            self.junction_table.setItem(i, 2, QCustomTableWidgetItemFloat(round(j.ppm, 3)))
            self.junction_table.setItem(i, 3, QtWidgets.QTableWidgetItem(j.orf.replace("_", " ").title()))
            self.junction_table.setItem(i, 4, QtWidgets.QTableWidgetItem(j.frame.replace("_", " ").title()))
        self.junction_table.setSortingEnabled(True)
        self.junction_table.sortItems(0, 0)

    @QtCore.pyqtSlot()
    def on_junction_table_selection_changed(self):
        selected = self.junction_table.selectionModel().selectedRows()
        if not selected:
            self.extrapolate_btn.setEnabled(False)
            self.selected_junction = None
            return
        row = selected[0].row()
        idx = self.junction_table.item(row, 0).data(QtCore.Qt.UserRole)
        direction, j = self.table_row_data[idx]
        self.selected_junction = (direction, j)
        self.rd_plot.mark_position(j.position)
        self.fivep_plot.mark_position(j.position)
        self.threep_plot.mark_position(j.position)
        self.extrapolate_btn.setEnabled(True)

    @QtCore.pyqtSlot()
    def on_rd_read_list_finished(self):
        self.calculate_depth_thread.start()

    @QtCore.pyqtSlot()
    def on_rd_depth_calculated(self):
        self.rd_plot.plot(self.results, self.selected_gene['orf_start'], self.selected_gene['orf_stop'])
        self.rd_load_btn.setEnabled(True)
        self.status_bar.showMessage("Read depth calculated.", 3000)

    # ---- Extrapolate Junction ----

    @QtCore.pyqtSlot()
    def on_extrapolate_btn_clicked(self):
        if not self.selected_junction:
            return
        direction, j = self.selected_junction
        base = self.fivep_base if direction == '5p' else self.threep_base
        file_prefix = self.fivep_file_prefix if direction == '5p' else self.threep_file_prefix
        if base is None:
            return
        blast_txt_path = os.path.join(self.directory, 'blast_results', file_prefix + base + '.blast.txt')
        junctions_txt_path = os.path.join(self.directory, 'junction_files', file_prefix + base + '.junctions.txt')
        if not os.path.exists(blast_txt_path) or not os.path.exists(junctions_txt_path):
            QtWidgets.QMessageBox.warning(self, "Extrapolate Junction",
                                          "Could not find the raw BLAST/junction files for this dataset:\n%s\n%s" %
                                          (blast_txt_path, junctions_txt_path))
            return

        self.extrapolate_btn.setEnabled(False)
        self.progress_dialog = QtWidgets.QProgressDialog("Scanning raw BLAST results for this junction...",
                                                          "Cancel", 0, 100, self)
        self.progress_dialog.setWindowTitle("Extrapolate Junction")
        self.progress_dialog.setWindowModality(QtCore.Qt.WindowModal)
        self.progress_dialog.setMinimumDuration(0)
        self.progress_dialog.setValue(0)

        self.extrapolate_thread = ExtrapolateJunctionThread(blast_txt_path, junctions_txt_path, direction,
                                                             self.selected_accession_number, j.position, j.query_start)
        self.extrapolate_thread.progress.connect(self.on_extrapolate_progress)
        self.extrapolate_thread.finished_result.connect(self.on_extrapolate_finished)
        self.progress_dialog.canceled.connect(self.extrapolate_thread.cancel)
        self.extrapolate_thread.start()
        self.progress_dialog.show()

    @QtCore.pyqtSlot(int, str)
    def on_extrapolate_progress(self, pct, message):
        self.progress_dialog.setValue(pct)
        self.progress_dialog.setLabelText(message)

    @QtCore.pyqtSlot(object, object, object, object, str)
    def on_extrapolate_finished(self, read_id, junction_seq, full_read, is_backwards, direction):
        self.progress_dialog.close()
        self.extrapolate_btn.setEnabled(True)
        if read_id is None:
            self.junction_seq_txt.setPlainText("No supporting read found for this junction.")
            return
        label = "5'" if direction == '5p' else "3'"
        strand = "Backward" if is_backwards else "Forward"
        fasta = ">%s %s %s\n%s\n\n>%s %s\n%s" % (self.selected_gene_name, label, strand, junction_seq,
                                                  self.selected_gene_name, read_id, full_read)
        self.junction_seq_txt.setPlainText(fasta)

    @QtCore.pyqtSlot()
    def on_copy_btn_clicked(self):
        QtWidgets.QApplication.clipboard().setText(self.junction_seq_txt.toPlainText())
        self.status_bar.showMessage("Copied to clipboard", 3000)


def appExit():
    app.quit()
    sys.exit()


if __name__ == '__main__':
    form = FragFinder_Gui(main_directory, gene_list_file, combined)
    form.show()
    app.aboutToQuit.connect(appExit)
    app.exec_()
