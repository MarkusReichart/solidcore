import collections, glob, json, os, sys

root = sys.argv[1]
def load(d):
    out = {}
    for p in glob.glob(os.path.join(root, d, '*.json')):
        u = os.path.basename(p)[:-5]
        r = json.load(open(p, encoding='utf-8'))
        out[u] = {c.get('name'): c for c in r.get('cases', [])}
    return out
on, off = load('reports'), load('reports-none')
common = set(on) & set(off)
print('units on', len(on), 'off', len(off), 'common', len(common))
tot = collections.Counter()
cat = collections.defaultdict(collections.Counter)
samples = {}
def fails(c):
    return c is not None and c.get('pass') is False
def ms(c):
    m = c.get('mismatches') or c.get('diff') or []
    return m if isinstance(m, list) else m.get('mismatches', [])
def rust_refuses(c):
    for m in ms(c):
        if 'status' in str(m.get('field', '')):
            a, b = str(m.get('cpp')), str(m.get('rust'))
            ok = ('0', '0e0', '1000', '1e3')
            return a in ok and b not in ok
    return False
for u in common:
    for name in set(on[u]) | set(off[u]):
        a, b = on[u].get(name), off[u].get(name)
        fa, fb = fails(a), fails(b)
        if fa and fb:
            k = 'fail_both'
        elif fa:
            k = 'flags_break_it'
            if rust_refuses(a):
                k = 'flags_break_it_rust_refuses'
        elif fb:
            k = 'flags_fix_it'
        else:
            k = 'pass_both'
        tot[k] += 1
        cat[u][k] += 1
        if k != 'pass_both':
            samples.setdefault((u, k), []).append(name)
print(tot.most_common())
for k in ['flags_break_it_rust_refuses', 'flags_break_it', 'fail_both', 'flags_fix_it']:
    print('==', k)
    rows = sorted(((v[k], u) for u, v in cat.items() if v[k]), reverse=True)[:24]
    for n, u in rows:
        print(n, u, [str(x)[:38] for x in samples[(u, k)][:2]])
json.dump({f'{u}|{k}': v for (u, k), v in samples.items()}, open(os.path.join(root, 'classes.json'), 'w'))
