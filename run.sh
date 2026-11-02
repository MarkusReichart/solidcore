set -u
pkill -f forkdiff/run10.sh; pkill -f 'target-forkdiff/difftest/difftest'; sleep 2
B=$HOME/solidcore-1808; E=/mnt/f/Dev/solidcore/evidence/2026-11-01-forkdiff; R=$E/r11; W=$B/wt-fd2; F=$B/fork-5e8c9c1
export CARGO_HOME=$B/cargo RUSTUP_HOME=$B/rustup PATH=$B/cargo/bin:$PATH CARGO_TARGET_DIR=$B/target-fd2 SOLIDCORE_CPP_DIR=$B/forkdiff-cpp-5e8c9c1
export RUSTFLAGS="-C relro-level=partial -C link-arg=-Wl,-z,lazy -C link-arg=-Wl,--unresolved-symbols=ignore-all"
mkdir -p $R/reports $R/reports-none
SHA=$(git -C /mnt/f/Dev/solidcore/repo rev-parse upstream/main); cd $B/repo; git fetch -q /mnt/f/Dev/solidcore/repo "$SHA" || exit 1
[ -d $W ] || git worktree add -q --detach $W "$SHA"; cd $W; git checkout -q -f --detach "$SHA"; git clean -qfd harness crates; git log -1 --format='rust %h %cd %s' --date=short | cut -c1-90; echo "fork $(git -C $F log -1 --format=%h)"
rm -rf upstream; mkdir upstream; ln -s $F/source upstream/source; ln -s $B/repo/upstream/TestFiles upstream/TestFiles; ln -s $B/repo/upstream/tests upstream/tests
python3 - <<'PY'
p='crates/solidcore/src/kernel.rs'
s=open(p,newline='\n').read()
a='fixes: UpstreamFixes::default(),'
if s.count(a)!=1:
    import re,sys
    print('kernel patch: pattern count',s.count(a)); print(re.findall(r'fixes: [^\n]*',s)[:6]); sys.exit(3)
b='''fixes: {
                let mut f = if std::env::var_os("FORKDIFF_NONE").is_some() { UpstreamFixes::default() } else { UpstreamFixes::all() };
                for (var, on) in [("FORKDIFF_OFF", false), ("FORKDIFF_ON", true)] {
                    if let Ok(v) = std::env::var(var) {
                        for n in v.split(',').filter(|n| !n.is_empty()) {
                            if let Some(x) = f.flag_mut(n) { *x = on; }
                        }
                    }
                }
                f
            },'''
open(p,'w',newline='\n').write(s.replace(a,b))
PY
[ $? = 0 ] || exit 1
: > $R/removed.txt
python3 - <<'PY'
import os,re
E='/mnt/f/Dev/solidcore/evidence/2026-11-01-forkdiff'
for n,dst in [('sync_copies.py','/tmp/sync11.py'),('restore_small.py','/tmp/restore11.py')]:
    s=open(f'{E}/{n}').read().replace('/wt-forkdiff','/wt-fd2')
    s=s.replace("B + '/fork-compare-2026-11-01/fork/source'","B + '/fork-5e8c9c1/source'")
    s=s.replace("R = B + '/repo/harness/cpp/shim'","R = B + '/wt-fd2/harness/cpp/shim'")
    open(dst,'w').write(s)
PY
python3 /tmp/sync11.py | head -1
for n in blend_fillet_corner_p3 blend_fillet_executive_p3 blend_fillet_geom_p1 blend_fillet_solver_p1 blend_fillet_solver_p1_static blend_fillet_solver_p2 blend_fillet_standard_solver_p1 primitive_primitive_creation_p1; do
  [ -f harness/cpp/shim/$n.cpp ] && { cp harness/cpp/shim/$n.cpp /tmp/removed_$n.cpp; rm -f harness/cpp/shim/$n.cpp; echo $n >> $R/removed.txt; }
done
python3 $E/stubs.py $E/missing.txt
python3 /tmp/restore11.py 2>&1 | tail -1 | cut -c1-110
for try in 1 2 3 4 5 6 7 8 9 10 11 12; do
  nice -n 8 cargo build --profile difftest -j 8 -p difftest --features cpp > $R/build.log 2>&1; rc=$?
  echo "build try $try rc=$rc"
  [ $rc = 0 ] && break
  bad=$(grep -E 'shim/[a-z0-9_]+\.cpp:[0-9]+:[0-9]+: (fatal )?error' $R/build.log | sed -E 's/.*shim\/([a-z0-9_]+)\.cpp.*/\1/' | sort -u)
  if [ -z "$bad" ]; then grep -E 'error' $R/build.log | sort -u | cut -c1-170 | head -12; exit 1; fi
  for n in $bad; do
    echo "shim out: $n ($(grep -c "shim/$n.cpp:[0-9]*:[0-9]*: error" $R/build.log) errors): $(grep -E "shim/$n\.cpp:[0-9]+:[0-9]+: (fatal )?error" $R/build.log | sed -E 's/^.*error: //' | cut -c1-90 | head -1)"
    echo $n >> $R/removed.txt; cp harness/cpp/shim/$n.cpp /tmp/removed_$n.cpp; rm -f harness/cpp/shim/$n.cpp
    grep -oE 'sc_[a-z0-9_]+' /tmp/removed_$n.cpp | sort -u >> $R/extra_syms.txt
  done
  # the removed shims' entry points become aborting stubs, unless another shim defines them
  python3 - <<'PY'
import glob,re,os
R='/mnt/f/Dev/solidcore/evidence/2026-11-01-forkdiff/r11/extra_syms.txt'
p='harness/cpp/shim/zz_forkdiff_stubs.cpp'
s=open(p).read()
have=set(re.findall(r'extern "C" void (\w+)\(\)',s))
others=''.join(open(f,errors='replace').read() for f in glob.glob('harness/cpp/shim/*.cpp') if not f.endswith('zz_forkdiff_stubs.cpp'))
ffi=''.join(open(f,errors='replace').read() for f in glob.glob('crates/solidcore-ffi/src/**/*.rs',recursive=True))
add=[]
for n in sorted(set(open(R).read().split())):
    if n in have or re.search(r'\b'+n+r'\s*\(',others) or not re.search(r'\b'+n+r'\b',ffi): continue
    add.append(n)
open(p,'a').write(''.join('extern "C" void %s() { stub("forkdiff stub %s"); }\n'%(n,n) for n in add))
print('stubs added',len(add))
PY
done
[ $rc = 0 ] || exit 1
BIN=$CARGO_TARGET_DIR/difftest/difftest
timeout 2400 nice -n 8 $BIN --green --jobs 8 --report-dir $R/reports > $R/green.log 2>&1; echo "on=$?"; grep -E '^difftest' $R/green.log | cut -c1-120
FORKDIFF_NONE=1 timeout 2400 nice -n 8 $BIN --green --jobs 8 --report-dir $R/reports-none > $R/green-none.log 2>&1; echo "none=$?"; grep -E '^difftest' $R/green-none.log | cut -c1-120
python3 $E/compare.py $R | sed -n 1,2p | cut -c1-220
