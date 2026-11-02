import sys
names = [l.strip() for l in open(sys.argv[1]) if l.strip()]
with open('harness/cpp/shim/zz_forkdiff_stubs.cpp', 'w') as f:
    f.write('#include <cstdio>\n#include <cstdlib>\n')
    f.write('static void stub(const char* n) { std::fputs(n, stderr); std::abort(); }\n')
    for n in names:
        f.write('extern "C" void %s() { stub("forkdiff stub %s"); }\n' % (n, n))
