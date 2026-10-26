# 4. Sort and sweep, rebuilt every step

Date: 2026-10-26

Status: Accepted

## Context

The broadphase needs to find pairs whose bounds overlap. The usual choice is a dynamic
AABB tree with fat bounds that update only when a body leaves them. That is fast for large
worlds, but the tree's shape depends on insertion order and history, and it is a lot of
code to keep deterministic and correct.

## Decision

Sort every shape's bounds on x each step and sweep. Bounds grow by the speculative distance
plus how far the body can travel in the step. A pair is kept only when at least one body is
dynamic and at least one is moving. Contacts between two bodies that are both still are
carried over untouched.

## Consequences

- Short, easy to verify, and deterministic.
- Sleeping piles cost a sort and a sweep, but no narrowphase.
- O(n log n) every step. Fine up to a few thousand shapes. A tree can replace it later
  without changing anything else.
