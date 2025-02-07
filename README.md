# Erdős–Straus Conjecture

This repository contains the code which I covered in a [blog post](https://www.anothermathblog.com/?p=1047) about the Erdős–Straus conjecture. You'll need to provide your own "primes.bin" to use that feature. I'm using a binary file containing the first 50 million primes stored as `std::uint32_t`, so a bit of modification on your part will be needed to support integers larger than `UINT32_MAX`.

Uses the [Number Theory Library](https://libntl.org/) by [Victor Shoup](https://www.shoup.net/).