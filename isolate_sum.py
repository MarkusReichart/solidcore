"""Isolation of flags_break_it: which flags, in how many units and samples; do their docs cite the fork?"""
import collections, json, re, subprocess
E = 'F:/Dev/solidcore/evidence/2026-11-01-forkdiff'
d = json.load(open(E + '/isolate-flags_break_it.json'))
kernel = subprocess.run(['git', '-C', 'F:/Dev/solidcore/repo', 'show', 'upstream/main:crates/solidcore/src/kernel.rs'],
                        capture_output=True, encoding='utf-8').stdout.split('\n')


def doc(flag):
    for i, l in enumerate(kernel):
        if re.match(r'\s*pub ' + flag + r': bool', l):
            j = i
            while kernel[j - 1].strip().startswith('///'):
                j -= 1
            return ' '.join(x.strip()[4:] for x in kernel[j:i])
    return ''


status = collections.Counter()
single = collections.defaultdict(set)
samples = collections.Counter()
combos = collections.Counter()
for k, v in d.items():
    u = k.split('|')[0]
    if isinstance(v, str):
        status[v[:30]] += 1
        continue
    status['isolated'] += 1
    key = '+'.join(sorted(v))
    combos[key] += 1
    single[key].add(u)
print('samples', len(d), status.most_common())
print('flag sets', len(combos))
for key, n in combos.most_common(45):
    fork = 'FORK' if all(re.search(r'fork|smlib@|tangency-labs/smlib', doc(f)) for f in key.split('+')) else '-'
    print(n, fork, key[:95], sorted(single[key])[:4])
