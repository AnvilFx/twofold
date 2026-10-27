# 10. Moved anchors wake what they overlap

Date: 2026-10-27

Status: Accepted

## Context

Waking follows contacts. A static or kinematic body that is moved wakes the dynamic bodies
touching it. But a static box teleported into a sleeping box it never touched woke
nothing. There was no contact to follow, and the broadphase skips pairs where neither body
is moving. After 60 steps the box was still asleep and still overlapping.

Options:

- Flag moved anchors, and let the next broadphase pair them with sleeping bodies. The flag
  has to survive until the next step, and the broadphase rules get a special case.
- Wake every sleeping body whose bounds overlap the anchor's new bounds right away, as
  Box2D does.

## Decision

When waking starts from a static or kinematic body, it also wakes every sleeping dynamic
body whose bounds, grown by the speculative distance, overlap that body's bounds and pass
the collision filters. A half-space is tested against its plane, since its bounds cover
all of space.

Removing a shape or body wakes only the dynamic bodies that touched it. Before, a floor
in that list woke everything else resting on it.

## Consequences

- `setTransform`, `setGeometry`, adding shapes and changing motion type are all covered by
  the same code.
- The cost is a pass over sleeping bodies when a static or kinematic body is changed
  directly. Kinematic bodies moved by velocity or targets don't pay it.
- Moving a ground plane wakes only what it reaches.
