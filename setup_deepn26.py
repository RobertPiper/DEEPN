r"""
py2app setup script for DEEPN_26v6 (Python 3 / PyQt5 / native Apple Silicon port).

Builds one app bundle containing:
  - deepn.py             (main launcher, Contents/MacOS/DEEPN_26v6)
  - gene_count_gui.py     \
  - junction_make_gui.py   \
  - gc_jm.py                \
  - mapster_gc_jm.py         \  extra_scripts: separate executables in
  - query_blast_gui.py        \  Contents/MacOS/, sharing the same embedded
  - read_depth_gui_v2.py       /  Python runtime as the main app.
  - fragfinder_gui.py         /
  - stat_maker_gui_v3.py     /

Stat Maker needs no state from DEEPN's own database (it does its own
folder/dataset selection internally), so it's bundled here as a launcher
module *and* still built separately as its own standalone app via
setup_statmaker26v3.py - same script, two independent builds, no shared
state between them. stat_maker_gui_v2.py (DESeq2-only, no Bayesian/MCMC
phase, no Bait2 specificity contrast) still exists in the repo and can
still be run directly, but is no longer what DEEPN's own Stat Maker
button launches, and is no longer bundled into this app.

mapster_gc_jm.py chains a bundled, self-contained copy of MAPster (built
separately - see LOCAL_MAPster/build_mapster.sh - and copied into
Contents/Resources/mapster/MAPster.app by build_deepn26.sh, since a nested
.app with its own Qt .framework symlink structure needs cp -R, not
py2app's data_files) into Gene Count then Junction Make, for the new
MAP+GC+JM button. MAPster's own genome selection and output folder are
pre-set from whatever library is selected in DEEPN before launch (see
_write_mapster_config() in deepn.py); everything else about running it -
picking input files, naming outputs, clicking Run - stays exactly as
manual as standalone MAPster already is.

Usage:
    python3 setup_deepn26.py py2app
"""
import os
import glob
from setuptools import setup

APP_NAME = 'DEEPN_26v6'

EXTRA_SCRIPTS = [
    'gene_count_gui.py',
    'junction_make_gui.py',
    'gc_jm.py',
    'mapster_gc_jm.py',
    'query_blast_gui.py',
    'read_depth_gui_v2.py',
    'fragfinder_gui.py',
    'stat_maker_gui_v3.py',
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
DATA_FILES += find_data_dir('r_scripts', 'r_scripts')
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
        'matplotlib', 'csv', 're', 'shutil', 'subprocess',
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
