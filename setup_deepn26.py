"""
py2app setup script for DEEPN_26v3 (Python 3 / PyQt5 / native Apple Silicon port).

Builds one app bundle containing:
  - deepn.py            (main launcher, Contents/MacOS/DEEPN_26v3)
  - gene_count_gui.py    \
  - junction_make_gui.py  \  extra_scripts: separate executables in
  - gc_jm.py               / Contents/MacOS/, sharing the same embedded
  - query_blast_gui.py    /  Python runtime as the main app.
  - read_depth_gui.py    /

Stat Maker is intentionally excluded - it was always a standalone app,
never launched from the DEEPN launcher itself.

Usage:
    python3 setup_deepn26.py py2app
"""
import os
import glob
from setuptools import setup

APP_NAME = 'DEEPN_26v4'

EXTRA_SCRIPTS = [
    'gene_count_gui.py',
    'junction_make_gui.py',
    'gc_jm.py',
    'query_blast_gui.py',
    'read_depth_gui_v2.py',
]

APP = [{
    'script': 'deepn.py',
    'extra_scripts': EXTRA_SCRIPTS,
    'plist': {
        'CFBundleName': APP_NAME,
        'CFBundleDisplayName': 'DEEPN',
        'CFBundleGetInfoString': 'DEEPN Sequencing',
        'CFBundleIdentifier': 'edu.uiowa.robertpiper.deepn26',
        'CFBundleShortVersionString': '26.1',
        'CFBundleVersion': '2601',
        'NSHumanReadableCopyright': '(c) 2016-2026 Venkatramanan Krishnamani, Robert C. Piper, Mark Stammnes',
        'NSHighResolutionCapable': True,
    },
}]

def find_data_files(sources, targets, patterns):
    ret = {}
    for i, source in enumerate(sources):
        target = targets[i]
        pattern = os.path.join(source, patterns[i])
        for filename in glob.glob(pattern):
            if os.path.isfile(filename):
                targetpath = os.path.join(target, os.path.relpath(filename, source))
                path = os.path.dirname(targetpath)
                ret.setdefault(path, []).append(filename)
    return sorted(ret.items())

def find_data_dir(root_dir, target_root):
    """Recursively collect every file under root_dir, preserving structure."""
    ret = {}
    for dirpath, _dirnames, filenames in os.walk(root_dir):
        rel = os.path.relpath(dirpath, root_dir)
        target = os.path.join(target_root, rel) if rel != '.' else target_root
        files = [os.path.join(dirpath, f) for f in filenames if not f.startswith('.')]
        if files:
            ret[target] = files
    return sorted(ret.items())

DATA_FILES = []
DATA_FILES += find_data_files(['ui', 'ui/Windows'], ['ui', 'ui/Windows'], ['*.ui', '*.ui'])
DATA_FILES += find_data_files(['functions'], ['functions'], ['*.py'])
DATA_FILES += find_data_dir('dictionaries', 'dictionaries')
DATA_FILES += find_data_dir('lists', 'lists')
DATA_FILES += find_data_dir('ncbi_blast', 'ncbi_blast')
DATA_FILES += find_data_dir('icon', 'icon')
DATA_FILES += [
    ('.', ['DEEPN_db.sqlite3', 'Y2Hreadme.txt', 'DragDropListView.py']),
]

OPTIONS = {
    'argv_emulation': False,
    'includes': [
        'sip', 'PyQt5', 'PyQt5.QtCore', 'PyQt5.QtGui', 'PyQt5.QtWidgets',
        'glob', 'pickle', 'time', 'sys', 'os', 'pydoc', 'json', 'numbers',
        'hashlib', 'decimal', 'threading', 'itertools', 'pyqtgraph',
        'joblib', 'sortedcontainers', 'xlsxwriter', 'numpy', 'pandas',
        'matplotlib',
    ],
    'packages': ['functions'],
    'iconfile': 'icon/Icon.icns',
    'excludes': ['PyQt4', 'tkinter'],
    'plist': APP[0]['plist'],
}

setup(
    app=APP,
    name=APP_NAME,
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],
    author='Venkatramanan Krishnamani, Robert C. Piper, Mark Stammnes',
    data_files=DATA_FILES,
)
