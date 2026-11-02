"""File-static functions the harness shims copied from C++, compared with the fork's: a stale copy
compares the fork against old upstream code, not against the fork."""
import glob, json, os, re

B = os.path.expanduser('~/solidcore-1808')
SHIM = B + '/repo/harness/cpp/shim'
FORK = B + '/fork-compare-2026-11-01/fork/source'
PIN = B + '/repo/upstream/source'
E = '/mnt/f/Dev/solidcore/evidence/2026-11-01-forkdiff'
HEAD = re.compile(r'^\s*(?:static\s+)?(?:inline\s+)?[A-Za-z_][\w:<>,\s\*&]*?\b((?:[A-Z]\w*::)?(?:sm|Sm|sm_)\w*)\s*\(', re.M)
STATIC = re.compile(r'^static\s+[\w:<>,\s\*&]+?\b(\w+)\s*\(', re.M)


def strip(s):
    s = re.sub(r'/\*.*?\*/', '', s, flags=re.S)
    s = re.sub(r'//[^\n]*', '', s)
    return re.sub(r'\s+', '', s)


def body(s, start):
    i = s.find('{', start)
    semi = s.find(';', start)
    if i < 0 or (0 <= semi < i):
        return None
    d = 0
    for k in range(i, len(s)):
        if s[k] == '{':
            d += 1
        elif s[k] == '}':
            d -= 1
            if d == 0:
                return s[start:k + 1]
    return None


def index(root):
    out = {}
    for p in glob.glob(root + '/**/*.cpp', recursive=True):
        try:
            s = open(p, encoding='utf-8', errors='replace').read()
        except OSError:
            continue
        for m in STATIC.finditer(s):
            b = body(s, m.start())
            if b:
                out.setdefault(m.group(1), []).append((os.path.relpath(p, root), strip(b)))
    return out


fork, pin = index(FORK), index(PIN)
rows = []
for p in sorted(glob.glob(SHIM + '/*.cpp')):
    s = open(p, encoding='utf-8', errors='replace').read()
    for m in STATIC.finditer(s):
        n = m.group(1)
        b = body(s, m.start())
        if not b or n not in pin:
            continue
        sb = strip(b)
        pin_same = any(sb == x for _, x in pin[n])
        if n not in fork:
            rows.append((os.path.basename(p), n, 'gone in fork', pin_same))
            continue
        fork_same = any(sb == x for _, x in fork[n])
        pin_fork_same = any(a == c for _, a in pin[n] for _, c in fork[n])
        if not fork_same:
            rows.append((os.path.basename(p), n, 'fork changed' if not pin_fork_same else 'shim edited', pin_same))
print('flagged', len(rows))
for r in rows:
    print(r)
json.dump(rows, open(E + '/stale_copies.json', 'w'), indent=0)
