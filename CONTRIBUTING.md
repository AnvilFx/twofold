# Contributing

Anyone can contribute. If you find a way to make Twofold faster, more stable, more
accurate or easier to follow, a pull request is welcome. So is a bug report with a scene
that shows the problem.

## Before you start

Small fixes can go straight to a pull request. For anything larger, such as a new
feature, a change to the solver or a new broadphase, open an issue first and describe
the idea. AnvilFx sets the overall direction of the engine, so it's worth agreeing on
the approach before you spend time on it.

Read the goals at the top of [docs/design.md](docs/design.md). A change has to keep them:

- **Deterministic.** The same calls in the same order give bit-identical results, and a
  copied world steps exactly like the original. Don't iterate in an order that depends on
  memory addresses or hashing, and don't add threads.
- **No dependencies.** The library uses only the C++17 standard library. The viewer may
  use its own.
- **Readable.** No clever tricks that make a step hard to follow. Speed matters, but not
  at that cost.

If you think one of these goals should change, open an issue and make the case.

The rest of [docs/design.md](docs/design.md) explains how the engine works, and
[docs/decisions](docs/decisions/) explains why. Check the decisions before changing
something that looks odd. It may have been tried another way already.

## Building and testing

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release
```

To run only the tests whose names contain some text, pass it to the test program:

```
build/twofold_tests rod
```

With Visual Studio or Xcode generators the program is in `build/Release/` instead.

The viewer, `twofold_viewer`, is the quickest way to see what a change does. Its demo
scenes are in `viewer/scenes.cpp`, and its Step panel shows step time, energy and the
deepest penetration.

## Reporting a bug

Include:

- the smallest scene that shows it, as code if you can
- the `WorldSettings` you use, if they aren't the defaults
- what you expected and what happened
- your platform and compiler

## What a pull request needs

- **A test.** A fix comes with a test that fails without it. A new feature comes with
  tests for what it does and for bad input. Tests live in `tests/test_*.cpp` and use the
  `TEST`, `CHECK` and `CHECK_NEAR` macros from `tests/test.h`. A new test file also goes in
  the list in `CMakeLists.txt`.
- **Numbers for any claim.** If a change is faster, give the step time before and after,
  the scene, the number of frames and your machine. If it is more stable, show what you
  measured and how, such as penetration depth or how often bodies come to rest. Every
  claim is checked before the change is merged, so say exactly how to repeat it.
- **No warnings.** The library builds with `/W4` on MSVC and `-Wall -Wextra` elsewhere.
- **Passing CI.** CI builds and tests on Windows, Linux with GCC and Clang, and macOS.
  Floating point results differ a little between platforms, so a test that passes on
  one machine can still fail on another. Give tests some tolerance.
- **Updated docs.** If behavior changes, update [docs/design.md](docs/design.md) and the
  README where they describe it. A choice between real alternatives gets a new record in
  [docs/decisions](docs/decisions/). Follow the format of the existing records, and
  replace an old record with a new one rather than editing it. A new limit you know of
  goes under Known limits in the design document.

Keep each pull request to one change. Describe what it changes and why.

The [devlog](docs/devlog.md) is kept by AnvilFx. You don't need to add to it.

## Code style

Match the code around your change. In short:

- 4 spaces, no tabs, LF line endings
- braces on the same line, and always braces, even around a one-line body
- `camelCase` for functions and variables, `PascalCase` for types, `kName` for constants
- lines up to about 120 columns
- in the library, `maxf` and `minf` from `src/internal.h` in place of `std::fmax` and
  `std::fmin`

Write a comment only to explain what the code doesn't already make clear. Keep it short
and state how the code behaves now. Leave out how it came to be that way: that belongs in
the commit message or a decision record.

## AI tools

Twofold itself is written by Claude, so contributions made with AI tools are welcome.
Say in the pull request if you used one. Either way, you are responsible for the change
and for every claim you make about it.

## License

By contributing, you agree that your work is released under the [MIT license](LICENSE).
