#!/usr/bin/env bash
# Builds DEEPN_26v4.app end to end: py2app, then two post-build fixes py2app
# doesn't handle correctly on its own, then codesigning and verification.
#
# Usage:
#   ./build_deepn26.sh
#
# Requires the project venv to be active (source .venv/bin/activate) and
# Homebrew's xz package installed (for the liblzma fix).
set -euo pipefail
cd "$(dirname "$0")"

APP="dist/DEEPN_26v4.app"

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
PYEOF

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
DB=$(ls "$APP/Contents/Resources/ncbi_blast/db/"*.ndb | head -1)
DB="${DB%.ndb}"
(cd "$APP/Contents/Resources" && ./ncbi_blast/bin/osx/arm64/blastn \
    -query "$TMP_QUERY" -db "$(python3 -c "import os; print(os.path.relpath('$DB', '.'))")" \
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
