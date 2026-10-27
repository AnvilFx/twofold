# 11. Impact contacts get a second pass

Date: 2026-10-27

Status: Accepted

## Context

A 0.8 m box falling at 6 m/s onto a small sphere resting on the floor pushed the sphere
70 mm into the floor. Gauss-Seidel visits each contact once per pass. The box and sphere
contact pushes the sphere down, and the sphere and floor contact only catches up on the
next substep. By then the sphere is already deep.

Options:

- More substeps. Fixes it, but every contact in the world pays.
- A second pass over every contact. Same cost problem, a bit cheaper.
- A second pass only where something is hitting hard.

## Decision

At the start of the step, every pair with a point approaching faster than the restitution
threshold marks its dynamic bodies as struck. All active contacts of struck bodies are
solved a second time right after the biased pass, and again after the relax pass.

## Consequences

- The sphere case drops from 70 mm to 47 mm.
- Scenes at rest pay nothing. Rain, where something is always landing, costs about 16%
  more.
- Pyramid, Tower and Dominoes rest exactly as before.
- Only the contacts touching a struck body get the second pass. An impact on top of a tall
  stack still takes a few substeps to reach the bottom.
