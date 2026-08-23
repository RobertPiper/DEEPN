"""
py2app setup script for StatMaker_26v6 (Python 3 / PyQt5 / native Apple Silicon,
DESeq2 poscounts normalization + Bayesian/MCMC + hit-criteria panel with a
reloadable general .csv).

Standalone app - matches the original architecture where Stat Maker was always
its own separate .app, never bundled inside DEEPN.app.

Usage:
    python3 setup_statmaker26v6.py py2app
"""
import os
import glob
from setuptools import setup

APP_NAME = 'StatMaker_26v6'

APP = [{
    'script': 'stat_maker_gui_v6.py',
    'plist': {
        'CFBundleName': APP_NAME,
        'CFBundleDisplayName': 'Stat Maker',
        'CFBundleGetInfoString': 'DEEPN Stat Maker',
        'CFBundleIdentifier': 'edu.uiowa.robertpiper.statmaker26v6',
        'CFBundleShortVersionString': '26.6',
        'CFBundleVersion': '2606',
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
    ret = {}
    for dirpath, _dirnames, filenames in os.walk(root_dir):
        rel = os.path.relpath(dirpath, root_dir)
        target = os.path.join(target_root, rel) if rel != '.' else target_root
        files = [os.path.join(dirpath, f) for f in filenames if not f.startswith('.')]
        if files:
            ret[target] = files
    return sorted(ret.items())

DATA_FILES = []
DATA_FILES += find_data_files(['ui'], ['ui'], ['Stat_Maker_v3.ui'])
DATA_FILES += find_data_files(['functions'], ['functions'], ['*.py'])
DATA_FILES += find_data_dir('r_scripts', 'r_scripts')
DATA_FILES += find_data_dir('icon', 'icon')
DATA_FILES += [
    ('.', ['DragDropListView.py']),
]

OPTIONS = {
    'argv_emulation': False,
    'includes': [
        'sip', 'PyQt5', 'PyQt5.QtCore', 'PyQt5.QtGui', 'PyQt5.QtWidgets',
        'glob', 'pickle', 'time', 'sys', 'os', 'pydoc', 'json', 'numbers',
        'hashlib', 'decimal', 'threading', 'itertools',
        'xlsxwriter', 'csv', 're', 'shutil', 'subprocess',
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
