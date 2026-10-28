# 26. Events follow the solve, and sensors test overlap

Date: 2026-10-28

Status: Accepted

## Context

Code that needs to know when shapes meet, for sounds, damage or triggers, had to call
`World::contacts` every step and diff the result. That list has a contact per facing for
mesh pairs and holds speculative contacts that never touch. There was also no way to have a
shape that detects overlap without colliding.

Options for when a pair counts as touching:

- When it has a contact with points, as in Box2D v3. Stable and cheap, but bodies that only
  pass within the speculative distance begin and end a touch.
- When a point overlaps. Contacts are springs, so a resting box hovers around zero overlap
  and would begin and end over and over.
- When the solver pushed on a point during the step, or a point overlaps.

Options for sensors:

- Sensor contacts in the contact list that the solver skips. Waking, islands, the mid-step
  refresh and continuous collision all walk contacts and would each need to skip them.
- A pass after the step that tests each sensor against what its bounds reach, kept apart from
  contacts, as Box2D v3.1 does.

## Decision

Contact events:

- A pair of shapes begins touching when the solver pushed on any of its points in some
  substep, or a point overlaps. It stops once nothing pushed and every point is further apart
  than the linear slop. A mesh pair is one pair however many contacts it has.
- Each contact keeps whether its pair was touching. The state is decided after the solve,
  and a pair whose contacts disappear at the start of the step ends there.
- A hit is the fastest point of a pair that the solver pushed on, when it approached faster
  than `WorldSettings::hitThreshold` at the start of the step. At most one per pair a step.
- Destroying a shape or body, or a joint that stops its bodies colliding, ends their touching
  pairs at once, with the ids from before the call. Filter and motion type changes are left to
  the next step.
- Events are always on. There are no per-shape flags.

Sensors:

- A sensor is a shape flag. It has no contacts and no mass, and continuous collision, waking
  and ray casts skip it.
- After continuous collision, each sensor queries the tree with its bounds and runs GJK
  against every shape it reaches that passes the collision filter, on any other body.
- A pair where neither body is moving keeps its last result, unless a call moved, reshaped,
  added or refiltered a shape since the last step.

Events from calls between steps are reported with the next step, before its own. Each step's
events are ordered by shape index.

## Consequences

- With no sensors, stepping is bit for bit what it was. The bench scenes time the same within
  noise.
- Speculative contacts never begin a touch. A pair that touches within a single substep and is
  then pulled apart by the mid-step refresh can be missed.
- Separations are measured at the start of the step, so an end can arrive a step after the
  shapes parted.
- A shape slot reused right after a destroy starts with no touches, since the old pairs ended
  before it was filled.
- A pair that still touches after its filter changed reports nothing.
- Sensors detect static, kinematic and sleeping bodies without waking them. A trigger zone in
  a still level costs a tree query a step. A half-space sensor visits every shape.
- Sensors test true 3D overlap, also between planar bodies. Two meshes, or a mesh and a
  half-space, are never detected. Joints don't hide shapes from sensors.
