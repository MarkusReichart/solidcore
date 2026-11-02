"""Known non-findings of the fork comparison, taken out of a report's failures.

Each rule names a difference in what the harness observes, not in a result of the kernel:
- brep_data_undef: the fork stores SmTopologyData in an SmTObjArray, so its unset ids hold the
  constructor's SM_UNDEF_ULONG where upstream's calloc, and the port, hold 0.
- executive_tangency: the fork's creator holds half a degree; the port takes it in the API paths
  (use_creator_tangency_tol), the harness builds some executives directly and reads 2.

A case counts as explained only if every one of its mismatches falls under a rule.
"""
import collections, glob, json, os, sys


def rule(unit, m):
    f, c, r = str(m.get('field')), str(m.get('cpp')), str(m.get('rust'))
    if unit == 'topology.brep_data' and f.startswith('ints.ints[') and c == '123456789' and r == '0':
        return 'brep_data_undef'
    if unit == 'offset.offset_executive_p2' and f == 'values.executive.geometry_parameters[2]' and c == '5e-1' and r == '2e0':
        return 'executive_tangency'
    return None


def main(root):
    total, explained = 0, collections.Counter()
    left = collections.Counter()
    partly = collections.Counter()
    for p in glob.glob(os.path.join(root, 'reports', '*.json')):
        u = os.path.basename(p)[:-5]
        for c in json.load(open(p, encoding='utf-8')).get('cases', []):
            if c.get('pass') is not False:
                continue
            total += 1
            ms = c.get('mismatches') or []
            if isinstance(ms, dict):
                ms = ms.get('mismatches', [])
            rs = [rule(u, m) for m in ms]
            if ms and all(rs):
                explained[rs[0]] += 1
            else:
                left[u] += 1
                if any(rs):
                    partly[u] += 1
    print('failed', total, 'explained', sum(explained.values()), dict(explained), 'left', sum(left.values()))
    print('partly explained, still counted', dict(partly))
    print('top left', left.most_common(15))


main(sys.argv[1])
