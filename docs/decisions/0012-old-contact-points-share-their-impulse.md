# 12. Old contact points share their impulse

Date: 2026-10-27

Status: Accepted

## Context

A new contact point took the impulse of the first old point within the speculative
distance. When a manifold changed shape between steps, for example two corners of a box
becoming three, two new points could continue the same old point and both take its full
impulse. The warm start then pushed that side of the contact twice as hard as the step
before. Small stacks never settled and kept rocking.

## Decision

Each new point continues the nearest old point within the speculative distance. Old points
continued by more than one new point split their impulse evenly between them. The same
matching is used after the narrowphase and after the mid-step refresh.

## Consequences

- The warm start never adds impulse that wasn't there the step before.
- Measured over 108 seeds, Hulls left 16 bodies awake at the end instead of 38.
- On its own it doesn't stop the rocking in three box towers. See
  [decision 13](0013-relax-pass-takes-friction-in-reverse.md).
