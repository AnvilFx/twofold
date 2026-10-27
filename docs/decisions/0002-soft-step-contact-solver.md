# 2. Soft step contact solver

Date: 2026-10-26

Status: Accepted

## Context

Stacks have to stay up, piles have to come to rest, and mixed mass ratios (a heavy body on
a light one) are common. Overlap will happen and has to be recovered from without adding
energy.

Options looked at:

- Sequential impulses with Baumgarte stabilization, as in Box2D 2.4. Simple, but pushout
  adds energy and stacks need many iterations.
- Sequential impulses with split impulses or an NGS position pass. Better, but the
  position pass costs as much again.
- XPBD. Stiff and stable, but friction and restitution are harder to get right.
- The soft step from Box2D v3: substeps with soft constraints, a relax pass, and
  restitution applied at the end.

## Decision

Use the soft step. Each step runs a fixed number of substeps (4 by default), each one:
integrate velocities, warm start, solve with bias, integrate positions, relax. Contacts are
springs with a stiffness in hertz and a damping ratio, capped at a quarter of the substep
rate. Contacts against static, kinematic or sleeping bodies are twice as stiff. Pushout is
capped by `maxPushoutSpeed`. Restitution runs once after the substeps.

Separation between substeps comes from each body's motion since the start of the step, so
collision runs once per step.

## Consequences

- Stacks hold with one iteration per substep.
- Overlap recovers at a bounded speed and the relax pass removes the pushout velocity, so
  deep overlap does not launch bodies.
- Each substep costs a full solve.
- Contact points are fixed for the whole step. Fast spinning bodies outrun them, which is
  what [decision 8](0008-collide-fast-pairs-again-mid-step.md) deals with.
