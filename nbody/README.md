# N-body problem
Mostly single-threaded (There are some papers showing possible parallel solutions using GPUs).
See the [description](https://benchmarksgame-team.pages.debian.net/benchmarksgame/description/nbody.html "n-body description")
and the Wikipedia [article](https://en.wikipedia.org/wiki/N-body "N-body simulation").

## C implemenations
All of these are in standard C99 (except #9, it uses C23 features)
and don't use compiler directives. Versions 1-3 and 6 don't use
array parameter declarations, so MSVC can compile them.

Some things that didn't work:
 - Manually flattening loops by precalculating body pairs as pairs
 of indices. This confused GCC's optimizer.
 - Replacing some branches with branchless code. While GCC was able
 to perform some new optimizations, the program became about 10% slower.
 
GCC command used:
```sh
gcc -O3 -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wvla -Werror -pedantic-errors -pipe -march=native -ffast-math -flto=auto -std=c99 -ffp-contract=fast -fmerge-all-constants -fgcse-sm -fgcse-las -fext-dce -fira-hoist-pressure -fselective-scheduling -fsel-sched-pipelining -fsel-sched-pipelining-outer-loops -fipa-reorder-for-locality -fipa-pta -ffinite-loops -fgraphite-identity -floop-nest-optimize -ftree-loop-im -ftree-loop-ivcanon -fvariable-expansion-in-unroller -fallow-store-data-races -fstdarg-opt -o nbody nbodyX.c -lm
```

MSVC (cl.exe) command used:
```cmd
cl /O2 /options:strict /W4 /utf-8 /validate-charset /MP /arch:AVX2 /fp:fast /jumptablerdata /GL /Gw /Fe:nbody /std:c17 nbodyX.c
```
