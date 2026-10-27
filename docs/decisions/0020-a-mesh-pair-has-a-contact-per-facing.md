# 20. A mesh pair has a contact per facing

Date: 2026-10-27

Status: Accepted

## Context

A contact has one normal and at most four points. A box resting in the crease between two
mesh faces touches both, along two normals. One manifold for the whole pair would have to
pick one of them, and the box would slide down the other face.

## Decision

A contact gains a part index next to its key. A mesh pair keeps one contact per group of
triangles that face the same way, sorted by key and then part. Every other pair still has
exactly one.

After the narrowphase, the old contacts of a mesh pair are reused only when all of them
still fit. Each new contact takes its impulses from the old contact whose normal is
closest. The mid-step refresh collides a mesh pair once and gives each contact the new
manifold closest to its normal.

## Consequences

- A box in a crease rests on two contacts and stays put.
- Code that walks contacts sees several with the same two shapes. `World::contacts` returns
  each one.
- Pairs without a mesh step exactly as before.
