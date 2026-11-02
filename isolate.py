"""For cases that pass against the fork with no flag and fail with all: the flag(s) responsible."""
import concurrent.futures as cf, json, os, re, subprocess, sys

B = os.path.expanduser('~/solidcore-1808')
BIN = B + '/target-forkdiff/difftest/difftest'
E = '/mnt/f/Dev/solidcore/evidence/2026-11-01-forkdiff'
CWD = B + '/wt-forkdiff'
src = open(CWD + '/crates/solidcore/src/kernel.rs', encoding='utf-8').read()
names = re.findall(r'^\s*pub ([a-z0-9_]+): bool,', src[src.index('upstream_fixes! {\n/// Opt-in'):], re.M)
print('flags', len(names), flush=True)
classes = json.load(open(E + '/classes.json'))
KEY = sys.argv[1]
PER_UNIT = int(sys.argv[2])

def passes(unit, case, on):
    env = dict(os.environ, FORKDIFF_NONE='1', FORKDIFF_ON=','.join(on), SOLIDCORE_CPP_DIR=B + '/forkdiff-cpp')
    try:
        p = subprocess.run([BIN, '--unit', unit, '--filter', case, '--no-report'], cwd=CWD, env=env,
                           capture_output=True, text=True, errors='replace', timeout=240)
    except subprocess.TimeoutExpired:
        return None
    return p.returncode == 0

def minimise(unit, case):
    on = list(names)
    if passes(unit, case, []) is not True:
        return 'not clean with none'
    if passes(unit, case, on) is not False:
        return 'not failing with all'
    n = 2
    while len(on) > 1:
        chunk = max(len(on) // n, 1)
        removed = False
        i = 0
        while i < len(on):
            rest = on[:i] + on[i + chunk:]
            if passes(unit, case, rest) is False:
                on = rest
                removed = True
            else:
                i += chunk
        if not removed:
            if chunk == 1:
                break
            n = min(n * 2, len(on))
    return on

jobs = []
for k, v in classes.items():
    u, c = k.split('|')
    if c == KEY:
        for name in v[:PER_UNIT]:
            jobs.append((u, name))
print('jobs', len(jobs), flush=True)
results = {}
with cf.ThreadPoolExecutor(10) as ex:
    futs = {ex.submit(minimise, u, n): (u, n) for u, n in jobs}
    for f in cf.as_completed(futs):
        u, n = futs[f]
        try:
            r = f.result()
        except Exception as e:
            r = 'error ' + str(e)[:80]
        results[f'{u}|{n}'] = r
        print(u, str(n)[:44], r, flush=True)
json.dump(results, open(f'{E}/isolate-{KEY}.json', 'w'), indent=0)
