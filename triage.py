import collections, glob, json, os, re, sys

rep = sys.argv[1]
ERR = re.compile(r'status|err', re.I)
by_cat = collections.Counter()
per_unit = collections.defaultdict(collections.Counter)
samples = collections.defaultdict(list)
keys = collections.Counter()
total = 0
for p in glob.glob(os.path.join(rep, '*.json')):
    r = json.load(open(p, encoding='utf-8'))
    u = os.path.basename(p)[:-5]
    for c in r.get('cases', []):
        if c.get('pass') is not False:
            continue
        total += 1
        for k in c:
            keys[k] += 1
        ms = c.get('mismatches') or c.get('diff') or []
        if isinstance(ms, dict):
            ms = ms.get('mismatches', [])
        cat = 'value'
        first = ms[0] if ms else {}
        fld = str(first.get('field', ''))
        cpp = str(first.get('cpp', ''))
        rust = str(first.get('rust', ''))
        outcome = str(c.get('outcome', ''))
        if not ms:
            cat = 'no_mismatch:' + outcome[:40]
        elif any('status' in str(m.get('field', '')) for m in ms):
            st = next(m for m in ms if 'status' in str(m.get('field', '')))
            a, b = str(st.get('cpp', '')), str(st.get('rust', ''))
            ok = lambda s: s in ('0', '0e0', '1000', '1e3')
            if ok(a) and not ok(b):
                cat = 'rust_refuses_fork_builds'
            elif ok(b) and not ok(a):
                cat = 'fork_refuses_rust_builds'
            else:
                cat = 'both_fail_other_status'
        elif fld.endswith('.len'):
            cat = 'different_length'
        else:
            try:
                e = float(first.get('rel_err') or 0)
                cat = 'value_small' if e < 1e-6 else 'value_large'
            except Exception:
                cat = 'value'
        by_cat[cat] += 1
        per_unit[u][cat] += 1
        if len(samples[(u, cat)]) < 2:
            samples[(u, cat)].append((c.get('name', '')[:42], fld[:30], cpp[:22], rust[:22]))
print('failed cases', total)
print(by_cat.most_common())
print('keys', keys.most_common(12))
for cat in ['rust_refuses_fork_builds', 'fork_refuses_rust_builds', 'both_fail_other_status']:
    print('==', cat)
    rows = sorted(((v[cat], u) for u, v in per_unit.items() if v[cat]), reverse=True)[:22]
    for n, u in rows:
        print(n, u, samples[(u, cat)][0])
print('== value_large top')
rows = sorted(((v['value_large'], u) for u, v in per_unit.items() if v['value_large']), reverse=True)[:14]
for n, u in rows:
    print(n, u, samples[(u, 'value_large')][0])
