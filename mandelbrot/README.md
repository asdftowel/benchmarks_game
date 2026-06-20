# Mandelbrot
Draw a b/w Mandelbrot set of specified size and save it in
the portable bitmap format. For a complete description see
[mandelbrot description (Benchmarks Game)](https://benchmarksgame-team.pages.debian.net/benchmarksgame/description/mandelbrot.html "Task description").

## C programs
These try to be portable and only use the C standard library.
They were tested with MSVC in standards conformance mode
(`/permissive-`) and Clang with `-std=c11`.

Versions with the "-prune" suffix use the shortcut from Julia #7.

TODO:
 - [x] It seems that it is possible to remove the char size
 dependency by not relying on bitwise operations and replacing
 constants with macros such as `CHAR_BIT`.
 - [ ] Test whether batching IO by storing longs instead of chars
 actually improves performance.
