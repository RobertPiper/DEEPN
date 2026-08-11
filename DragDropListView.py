from PyQt5 import QtCore, QtGui, QtWidgets


class DragDropListView(QtWidgets.QListWidget):
    dropped = QtCore.pyqtSignal(str)
    deleted = QtCore.pyqtSignal(str)

    def __init__(self, type, parent=None):
        super(DragDropListView, self).__init__(parent)
        self.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.setIconSize(QtCore.QSize(72, 72))

    def startDrag(self, supportedActions):
        drag = QtGui.QDrag(self)
        t = self.selectedItems()[0].GetData()
        mimeData = self.model().mimeData(self.selectedIndexes())
        mimeData.setText(str(t))

        drag.setMimeData(mimeData)
        if drag.exec_(QtCore.Qt.MoveAction) == QtCore.Qt.MoveAction:
            for item in self.selectedItems():
                self.takeItem(self.row(item))

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.ignore()
        else:
            event.accept()

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls():
            event.ignore()
        else:
            event.accept()

    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
            text = event.mimeData().text()
            self.dropped.emit(text)
        else:
            event.ignore()

    def keyPressEvent(self, event):
        if event.key() == 16777219:
            self._del_item()

    def _del_item(self):
        for item in self.selectedItems():
            self.deleted.emit(item.GetData())
            self.takeItem(self.row(item))
