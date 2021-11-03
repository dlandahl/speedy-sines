Comparing the speed of a few possible ways to make sine oscillators.

An important disclaimer is that this is a contrived test due to frequency modulation not being implemented. In reality, you would just memcpy whole chunks of the signal if you didn't require real-time frequency changes. The performance of different techniques should still be mostly representative.

I use additive synthesis as a benchmark because it is a very powerful kind of synthesis that happens to require a large amount of sine function evaluations.

Most compilers automatically use fast approximations for libc sin when you clearly keep theta a small value, i.e. [-pi, pi] or [0, tau]. This means that -ffast-math makes no speed difference as long as the compiler can see the phase modulo happening before evaluating sin.

Lookup tables are still much better regardless. In my benchmark the speed is almost double. The compiler will not be able to vectorise them very well, however. This is not surprising because vpgatherdps and friends are quite specialised, as there is no special memory access hardware to go along with those instructions. If you touch different cache lines to fill a SIMD register, scatter / gather will not help that issue. But we can use them to make an oscillator, which is actually the fastest that I present here by a large margin. It loads and linearly interpolates eight float-32s at a time from a lookup table.

Since mathematical approximations can be vectorised more easily, if you don't have AVX2 specifically, then something like the AVX_MATHFUN library is considerably better than using a table. They are similar to the vector routines in SVML, which is the Intel-compiler exclusive vector maths library. Note however that the SSE based sine approximations from SSE_MATHFUN are a bit slower than a scalar lookup table.

I tried a sine approximation from stackoverflow. This is not faster than a scalar lookup table, and these kinds of routines tend to be bad for audio. It does not sound perfect and you can find spectral artifacting in a spectrum analyser.

Notice that I'm using a 64bit phase value. This is needed to prevent phase drift between different partials. Keeping phase values separately in each SIMD channel also causes too much drift.

Using the custom hand-written vpgatherdps-based AVX2 lookup table gives an order of magnitude speed improvement over the naive -ffast-math libc sine and a 5x speedup over a scalar lookup table.
