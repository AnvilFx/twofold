# 36. Speeds are limited and joints keep up

Date: 2026-11-02

Status: Accepted. Amends [decision 31](0031-world-inertia-is-fixed-for-the-step.md).

## Context

A fuzzer that creates, edits and steps random worlds reached NaN in 226 of 300 seeds with
locks and joints on, and in 19 of 300 with both off. Three things fed it.

A body of 1e-5 kg hung by a hinge from a 100 kg body, struck for half a second, took
impulses against the heavy body that grew every substep until they overflowed. Joint
impulses warm start, so a step that ended fast started the next one fast.

A distance joint holding a planar body to a point outside its plane can't be satisfied.
The locks leave the body almost no motion along the joint, so its effective mass there is
huge, and each correction threw the body past the closest distance it could reach. Trapped
bodies reached 700 m/s and 2000 rad/s.

A long body spinning fast turns a lot within a step, and its world inertia was fixed at
the start of the step. Contacts rebuilt at a mid-step refresh took lever arms from the
current pose and inertia from the old one. Spinning rods gained up to 77% of their energy.

Options:

- Clamp the speeds and leave the joints alone. That made joint fuzzing worse, 256 seeds
  against 193 in the tamer mode: a clamped body keeps its stored joint impulses, which
  drive it past the limit again from the next warm start and grow.
- Clear a joint's impulses whenever a clamp touches one of its bodies. Every clamp then
  drops the joint for a substep.
- Clamp the speeds and scale each joint's stored impulses by the fraction of speed its
  bodies kept.

## Decision

`WorldSettings` has `maxLinearSpeed` and `maxAngularSpeed`, both 1000. Dynamic body
velocities are scaled down to them after gravity and forces, before positions move, and
after restitution. Each body keeps the smallest fraction of its speed it kept since its
joints last warm started, and before each warm start the joints scale their stored impulses
by the smaller fraction of their two bodies, then reset it.

A distance joint scales its error by the share of the bodies' unlocked response along the
line that the locks leave. A limit that is not violated keeps its slack exact, and the
settle pass scales its moves the same way.

At a mid-step refresh every awake dynamic body refreshes its world inverse inertia before
the contacts are rebuilt. That inertia is still fixed between refreshes, and joint rows
build their axes as they solve, so they read the same inertia as the contacts.

## Consequences

- NaN seeds out of 300: 226 to 5 with locks and joints, 193 to 0 in the tamer mode, 214 to 2
  with locks only, 19 to 0 with neither. Joints only reach 1.
- The struck hinge settles where it went to infinity before. Trapped distance joints settle
  at about 1 m/s.
- Spinning rods gain at most 8% of their energy from 150 to 3000 rad/s.
- Anything that should move faster than the limits loses momentum and energy. `maxLinearSpeed`
  scales with `lengthUnit`, and a scene that needs more raises the limits.
- A planar body whose distance joint is exactly as long as its offset from the plane
  wobbles about 6 cm instead of being held stiffly.
