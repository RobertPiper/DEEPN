#!/usr/bin/env bash
# Builds DEEPN_26v6.app end to end: py2app, then two post-build fixes py2app
# doesn't handle correctly on its own, then codesigning and verification.
#
# Usage:
#   ./build_deepn26.sh
#
# Requires the project venv to be active (source .venv/bin/activate) and
# Homebrew's xz package installed (for the liblzma fix).
set -euo pipefail
cd "$(dirname "$0")"

APP="dist/DEEPN_26v6.app"

echo ">>> Cleaning previous build..."
rm -rf "$APP"

echo ">>> Running py2app..."
python3 setup_deepn26.py py2app

echo ">>> Fixing liblzma.5.dylib (py2app's copy has a stale/corrupted code"
echo "    signature that fails strict codesign validation)..."
LIBLZMA_SRC=$(find /opt/homebrew/Cellar/xz -name 'liblzma.5.dylib' | sort -V | tail -1)
if [ -z "$LIBLZMA_SRC" ]; then
    echo "ERROR: could not find a Homebrew liblzma.5.dylib under /opt/homebrew/Cellar/xz" >&2
    exit 1
fi
cp "$LIBLZMA_SRC" "$APP/Contents/Frameworks/liblzma.5.dylib"
install_name_tool -id "@executable_path/../Frameworks/liblzma.5.dylib" "$APP/Contents/Frameworks/liblzma.5.dylib"
codesign --force --sign - "$APP/Contents/Frameworks/liblzma.5.dylib"

echo ">>> Fixing ncbi_blast binaries' dylib paths..."
echo "    py2app rewrites every bundled binary's dylib references to"
echo "    @executable_path/../Frameworks/ and copies the actual dylibs to"
echo "    Contents/Frameworks/ - correct only for executables directly in"
echo "    Contents/MacOS/. blastn/blastdbcmd/makeblastdb live 4 directories"
echo "    deeper (Contents/Resources/ncbi_blast/bin/osx/arm64/), so that"
echo "    relative path resolves to the wrong place and every BLAST search"
echo "    silently fails to launch. Fix: copy the full transitive dylib"
echo "    closure into the location those binaries actually look in"
echo "    (ncbi_blast/bin/osx/Frameworks/), where @executable_path/../Frameworks/"
echo "    from arm64/ correctly resolves."
python3 - "$APP" <<'PYEOF'
import subprocess, os, shutil, sys

app = sys.argv[1]
bin_dir = os.path.join(app, 'Contents/Resources/ncbi_blast/bin/osx/arm64')
src = os.path.join(app, 'Contents/Frameworks')
dst = os.path.join(app, 'Contents/Resources/ncbi_blast/bin/osx/Frameworks')
os.makedirs(dst, exist_ok=True)

def deps_of(path):
    out = subprocess.check_output(['otool', '-L', path]).decode()
    result = set()
    for line in out.splitlines()[1:]:
        line = line.strip()
        if line.startswith('@executable_path/../Frameworks/'):
            result.add(line.split()[0].replace('@executable_path/../Frameworks/', ''))
    return result

queue = []
for binname in ('blastn', 'blastdbcmd', 'makeblastdb'):
    queue.extend(deps_of(os.path.join(bin_dir, binname)))

seen = set()
missing = []
while queue:
    dylib = queue.pop()
    if dylib in seen:
        continue
    seen.add(dylib)
    src_path = os.path.join(src, dylib)
    if not os.path.exists(src_path):
        missing.append(dylib)
        continue
    shutil.copy2(src_path, os.path.join(dst, dylib))
    queue.extend(d for d in deps_of(src_path) if d not in seen)

if missing:
    print("ERROR: missing source dylibs: %s" % missing, file=sys.stderr)
    sys.exit(1)
print("    Resolved and copied %d dylibs (transitive closure)." % len(seen))

# Some dylibs (seen in the wild: mbedtls's libmbedx509/libmbedcrypto) reference
# their own dependencies via bare @rpath/<soname> instead of
# @executable_path/../Frameworks/<soname>, with an LC_RPATH that points
# nowhere in our bundle (e.g. @loader_path/../lib, a directory that doesn't
# exist here). dyld then falls back to the absolute Homebrew path baked in
# at build time - which only "works" on a dev machine that happens to have
# the same Homebrew package installed, and hard-crashes (SIGABRT, "Library
# not loaded") on any other machine.
#
# Fix by redirecting each such dylib's LC_RPATH to '@loader_path' (its own
# directory, where the closure-copy above already placed every dependency)
# and adding a symlink under the exact soname each @rpath reference expects.
# Deliberately NOT using install_name_tool -change on the @rpath/<soname>
# dependency strings themselves: growing a load-command path string in
# place (the replacement is always longer) can silently corrupt the Mach-O
# header - confirmed by testing, it produced a dylib that made dyld hard-
# kill (SIGKILL, no diagnostic) partway through loading. Shrinking an
# LC_RPATH (the replacement here is always shorter) does not have this
# failure mode.
import re

def base_name(name):
    m = re.match(r'^(.*?)\.\d', name)
    return m.group(1) if m else name

def rpaths_of(path):
    out = subprocess.check_output(['otool', '-l', path]).decode().splitlines()
    paths = []
    for i, line in enumerate(out):
        if line.strip() == 'cmd LC_RPATH':
            m = re.search(r'path (\S+)', out[i + 2])
            if m:
                paths.append(m.group(1))
    return paths

fixed = 0
for d in sorted(seen):
    path = os.path.join(dst, d)
    out = subprocess.check_output(['otool', '-L', path]).decode()
    rpath_deps = [l.strip().split()[0][len('@rpath/'):]
                  for l in out.splitlines()[1:] if l.strip().startswith('@rpath/')]
    if not rpath_deps:
        continue
    for old_rpath in rpaths_of(path):
        subprocess.check_call(['install_name_tool', '-rpath', old_rpath, '@loader_path', path])
    # install_name_tool invalidates the embedded signature, and codesign
    # --deep on the whole app (run later in this script) does NOT reach
    # into this non-standard, 4-levels-deep Resources path - confirmed by
    # testing, it leaves the pre-modification signature in place, which
    # then legitimately mismatches the file content. At runtime the kernel
    # catches that mismatch and SIGKILLs the process (CODESIGNING, "Invalid
    # Page") the moment dyld tries to load it. Re-sign explicitly here,
    # same as the liblzma.5.dylib fix above does for the same reason.
    subprocess.check_call(['codesign', '--force', '--sign', '-', path])
    for ref_name in rpath_deps:
        base = base_name(ref_name)
        candidates = [f for f in seen if base_name(f) == base]
        if len(candidates) != 1:
            print("ERROR: cannot resolve @rpath/%s referenced by %s (candidates: %s)" % (ref_name, d, candidates), file=sys.stderr)
            sys.exit(1)
        link_path = os.path.join(dst, ref_name)
        if not os.path.exists(link_path):
            os.symlink(candidates[0], link_path)
    fixed += 1

print("    Fixed %d dylib(s) with unresolvable @rpath dependencies (redirected rpath, added soname symlinks)." % fixed)
PYEOF

echo ">>> Bundling MAPster..."
echo "    A nested .app with its own Qt .framework bundles has Versions/Current"
echo "    style symlinks that a naive per-file copy (py2app's data_files, or a"
echo "    plain Python file-walk) would break - cp -R preserves them correctly,"
echo "    so this is a straight directory copy, not routed through py2app."
MAPSTER_SRC="/Users/robertpiper2/LOCAL_MAPster/MAPster.app"
if [ ! -d "$MAPSTER_SRC" ]; then
    echo "ERROR: $MAPSTER_SRC not found - run LOCAL_MAPster/build_mapster.sh first" >&2
    exit 1
fi
mkdir -p "$APP/Contents/Resources/mapster"
rm -rf "$APP/Contents/Resources/mapster/MAPster.app"
cp -R "$MAPSTER_SRC" "$APP/Contents/Resources/mapster/MAPster.app"

echo ">>> Codesigning..."
codesign --deep --force --sign - "$APP"

echo ">>> Verifying..."
codesign --verify --deep --strict "$APP"

echo ">>> Smoke-testing blastn from inside the bundle..."
TMP_QUERY=$(mktemp -t deepn_build_test.XXXXXX.fa)
cat > "$TMP_QUERY" <<'EOF'
>test
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF
DB_NAME=$(basename "$(ls "$APP/Contents/Resources/ncbi_blast/db/"*.ndb | head -1)" .ndb)
(cd "$APP/Contents/Resources" && ./ncbi_blast/bin/osx/arm64/blastn \
    -query "$TMP_QUERY" -db "ncbi_blast/db/$DB_NAME" \
    -task blastn -dust no -num_threads 2 -outfmt 7 -out /tmp/deepn_build_smoketest.txt \
    -evalue 0.2 -max_target_seqs 10)
if [ -s /tmp/deepn_build_smoketest.txt ]; then
    echo "    blastn OK - produced output."
else
    echo "ERROR: blastn smoke test produced no output" >&2
    exit 1
fi
rm -f "$TMP_QUERY" /tmp/deepn_build_smoketest.txt

echo ""
echo ">>> Build complete: $APP"
