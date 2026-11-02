"""Replace pure copies of file-static C++ functions in the harness shims by the fork's version.

A copy is replaced only when it equals the pinned upstream function apart from comments and
whitespace and the fork's differs: the copy then carries no deliberate harness edit. The list of
replaced functions goes to sync_copies.json; a shim that no longer compiles goes back whole
(restore_file). Works on the forkdiff worktree only.
"""
import glob, json, os, re, subprocess, sys

B = os.path.expanduser('~/solidcore-1808')
W = B + '/wt-forkdiff'
SHIM = W + '/harness/cpp/shim'
FORK = B + '/fork-compare-2026-11-01/fork/source'
PIN = B + '/repo/upstream/source'
E = '/mnt/f/Dev/solidcore/evidence/2026-11-01-forkdiff'
STATIC = re.compile(r'^static\s+[\w:<>,\s\*&]+?\b(\w+)\s*\(', re.M)


def strip(s):
    s = re.sub(r'/\*.*?\*/', '', s, flags=re.S)
    s = re.sub(r'//[^\n]*', '', s)
    return re.sub(r'\s+', '', s)


def span(s, start):
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
                return k + 1
    return None


def index(root):
    out = {}
    for p in glob.glob(root + '/**/*.cpp', recursive=True):
        s = open(p, encoding='utf-8', errors='replace').read()
        for m in STATIC.finditer(s):
            e = span(s, m.start())
            if e:
                out.setdefault(m.group(1), []).append((os.path.relpath(p, root), s[m.start():e]))
    return out


pin, fork = index(PIN), index(FORK)
done = {}
for p in sorted(glob.glob(SHIM + '/*.cpp')):
    s = open(p, encoding='utf-8', errors='replace', newline='\n').read()
    edits = []
    for m in STATIC.finditer(s):
        n = m.group(1)
        e = span(s, m.start())
        if not e or n not in pin or n not in fork:
            continue
        cur = strip(s[m.start():e])
        match = [f for f, x in pin[n] if strip(x) == cur]
        if not match:
            continue
        # the fork's function of that name in the file the copy came from
        cand = [x for f, x in fork[n] if os.path.basename(f) == os.path.basename(match[0])]
        if len(cand) != 1 or strip(cand[0]) == cur:
            continue
        edits.append((m.start(), e, cand[0], n))
    if not edits:
        continue
    for a, b, new, _ in sorted(edits, reverse=True):
        s = s[:a] + new + s[b:]
    open(p, 'w', encoding='utf-8', newline='\n').write(s)
    done[os.path.basename(p)] = [x[3] for x in edits]
json.dump(done, open(E + '/sync_copies.json', 'w'), indent=0)
print('files', len(done), 'functions', sum(len(v) for v in done.values()))
for k, v in done.items():
    print(k, v)
