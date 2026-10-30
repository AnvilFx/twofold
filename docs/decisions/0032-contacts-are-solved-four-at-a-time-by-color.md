# 32. Contacts are solved four at a time by color

Date: 2026-10-30

Status: Accepted.

## Context

The contact solve took 49% of a step in Pyramid and 12% to 33% in the other scenes. Solved
one contact at a time, it waits on a chain of dependent float operations at each point: two
turns of the anchors, a square root, and a velocity read right after the last point wrote
it. Copying the bodies into a compact array first changed nothing, so memory isn't what
holds it back. Building with AVX2 changed nothing either, because the compiler doesn't
vectorize this code by itself.

Points of one contact share both bodies, so they can't be solved side by side. Contacts
that move different bodies can.

Options:

- Keep one contact at a time.
- Color the contacts as Box2D v3 does, so that no two contacts of one color move the same
  body, and solve each color several contacts at a time with SIMD.

## Decision

Each step, every active contact takes the lowest of 32 colors that neither of its moving
bodies already has. Static, kinematic and sleeping bodies take no colors. Within a color,
contacts with more points come first, and they are cut into batches of four with SSE2 or
NEON. Colors are solved in order, then planar contacts, contacts with rolling resistance
and contacts that found no free color, one at a time. The second pass for impact contacts
uses the same colors.

A batch keeps its contacts' normals, anchors and axes side by side, and builds the axes
itself at the start of the step and after every mid-step refresh. Warm starting and
restitution run on the batches too, so the per contact axes are only built for contacts
solved one at a time.

Each lane runs the float operations of the single contact solve in the same order. There
is no fused multiply-add and no approximate reciprocal or square root, maxima follow
`maxf` for NaN and signed zeros, and lanes or points a batch doesn't use are masked out
with selects. Building with `TWOFOLD_SIMD=OFF` solves the same batches one lane at a
time, and gives the same bits.

## Consequences

- Contacts are solved in color order instead of contact order, so every scene with
  contacts changes. Over 12 sets of 32 seeds against the previous commit:
  - Pyramid rests at 2.55 mm of overlap instead of 2.24 mm, and Tower at 7.49 mm instead
    of 7.29 mm. Contact order used to solve a stack roughly from the floor up.
  - Dominoes has more frames past 5 mm, 8995 against 6056, with the same worst overlap.
  - Rain has more frames past 50 mm and 100 mm, 1700 and 118 against 1478 and 66, and a
    lower resting maximum, 21 mm against 37 mm.
  - Meshes has a deeper worst overlap, 383 mm against 314 mm, with fewer frames past
    20 mm, 1750 against 1940, and a lower resting maximum, 18 mm against 24 mm. Its worst
    sets swing that much from small changes in the start of a GJK search.
  - Hulls, Sweeper, Joints and Sleeping are within noise. Cutouts, Flat layers and Slopes
    are unchanged, since their contacts are planar or few.
- Whole steps with sleep off, 12 seeds of 600 frames on one core, best of three: Pyramid
  0.480 to 0.329 ms, Hulls 0.585 to 0.542 ms, Meshes 1.705 to 1.601 ms, Tower 0.026 to
  0.021 ms, Sweeper 0.359 to 0.346 ms, Dominoes 0.133 to 0.128 ms. Rain, Joints and
  Cutouts are within 2%, since collision or joints take most of their step. Jolt runs
  Pyramid in 0.269 ms.
- Without SSE2 or NEON, for example WebAssembly without SIMD, the same batches run one lane
  at a time, about 1.7 times slower in the solve than the old single contact solve.
- A batch holds a copy of its contacts' setup, and the coloring and batching run every
  step.
