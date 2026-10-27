# 13. The relax pass takes friction in reverse

Date: 2026-10-27

Status: Accepted

## Context

Three box towers rocked and twisted at about 9 Hz and never went to sleep. The motion was
small but steady. Gauss-Seidel always solved a contact's friction points in the same
order, so the first point took up the error and the last one overshot, the same way every
substep. Friction kept feeding the motion.

Each of these stopped it:

- Solving friction only in the relax pass. Changes friction everywhere to fix one case.
- Removing the doubled stiffness on static contacts. Halves the stiffness against the
  ground, so every stack sags twice as far into it.
- Reversing the whole relax pass. The overlap tail in Rain got worse: 62 frames past
  100 mm instead of 29, over 108 seeds.
- Alternating the order every substep. Similar to the decision below, with a slightly
  worse tail.
- Reversing only the friction points in the relax pass.

## Decision

The biased pass solves friction points in order, the relax pass in reverse. Normal points
keep their order in both.

## Consequences

Over 108 seeds:

| Scene | Before | After |
| --- | --- | --- |
| Rain, frames past 50 / 100 / 150 mm | 527 / 32 / 5 | 441 / 29 / 0 |
| Flat layers, awake after 3600 frames | 16 | 5 |

- Towers, Cutouts and Hulls all end asleep.
- Tower, Pyramid, Dominoes and Slopes rest exactly as before.
- Costs nothing.
