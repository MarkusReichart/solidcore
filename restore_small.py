"""Bring six removed shims back into the forkdiff worktree, adapted to the fork's signatures."""
import os, re, shutil

B = os.path.expanduser('~/solidcore-1808')
W = B + '/wt-forkdiff/harness/cpp/shim'
R = B + '/repo/harness/cpp/shim'
names = ['offset_offset_geometry_creation', 'sweep_api_brep',
         'topology_brep_p10', 'geom_attribute', 'geom_a_object']


def sub(s, a, b, n=1):
    assert s.count(a) == n, a
    return s.replace(a, b)


for n in names:
    s = open(f'{R}/{n}.cpp', encoding='utf-8', newline='\n').read()
    if n == 'offset_offset_geometry_creation':
        s = sub(s, 'SmSSIData ssi;\n            SmStatus s = (op == S_LOCAL_FACEUSE_OFFSET)',
                'SmSSIData ssi;\n            SmBoolean collapses = FALSE; // forkdiff: the fork\'s out-parameter\n            SmStatus s = (op == S_LOCAL_FACEUSE_OFFSET)')
        s = sub(s, 'local_now.FaceuseOffset(ctx, fu, skip ? TRUE : FALSE, ns, ssi)',
                'local_now.FaceuseOffset(ctx, fu, skip ? TRUE : FALSE, ns, ssi, collapses)')
        s = sub(s, 'geom.FaceuseOffset(ctx, fu, skip ? TRUE : FALSE, ns, ssi)',
                'geom.FaceuseOffset(ctx, fu, skip ? TRUE : FALSE, ns, ssi, collapses)')
    if n == 'sweep_api_brep':
        s = sub(s, 'SmApiCreateTaperExtrude(curves, ha[0], ha[1], caps, corner ? TRUE : FALSE, res)',
                'SmApiCreateTaperExtrude(curves, ha[0], ha[1], caps, corner ? SM_OC_FILLET_CORNER : SM_OC_LINEAR_EXTENSION, res)')
    if n == 'topology_brep_p10':
        s = sub(s, 'd.m_vShells = shells;',
                'd.m_vShells.ReSet(); for (ULONG k = 0; k < shells.GetSize(); k++) d.m_vShells.Add(shells[k]);')
        s = sub(s, 'd.m_vRegions = regions;',
                'd.m_vRegions.ReSet(); for (ULONG k = 0; k < regions.GetSize(); k++) d.m_vRegions.Add(regions[k]);')
    if n in ('geom_attribute', 'geom_a_object'):
        s = sub(s, '#include "sc_shim.h"\n', '#include "sc_shim.h"\n#include <SmShell.h>\n')
    open(f'{W}/{n}.cpp', 'w', encoding='utf-8', newline='\n').write(s)

p = f'{W}/zz_forkdiff_stubs.cpp'
stubs = open(p, encoding='utf-8', newline='\n').read().split('\n')
text = ''.join(open(f'{W}/{n}.cpp', encoding='utf-8', errors='replace').read() for n in names)
keep = []
dropped = 0
for line in stubs:
    m = re.match(r'extern "C" void (\w+)\(\)', line)
    if m and re.search(r'\b' + m.group(1) + r'\s*\(', text):
        dropped += 1
        continue
    keep.append(line)
open(p, 'w', encoding='utf-8', newline='\n').write('\n'.join(keep))
print('restored', len(names), 'stubs dropped', dropped, 'stubs left', sum(l.startswith('extern') for l in keep))
