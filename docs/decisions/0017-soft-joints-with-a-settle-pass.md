# 17. Soft joints with a settle pass

Date: 2026-10-27

Status: Accepted

## Context

Joints have to fit the soft step the contacts already use. Solving them as springs in the
same substeps is simple and stable, but a spring stretches under load. A chain of 12 light
links holding a ball about 200 times heavier sags noticeably, and the error builds up
along the chain.

Options:

- Rigid joints with their own position correction, as in Box2D v2. A second solver
  structure beside the contacts.
- Soft joints only, as in Box2D v3. Simple, but the heavy chain stretches.
- Softer joints, so the springs fight the other rows less. At 30 Hz the still chain ends
  0.67 m long in total against 0.17 m at 60 Hz.
- Soft joints, then a short position pass that moves bodies to close what error is left,
  as in Jolt's position steps.

## Decision

Joints are solved in the substeps before the contacts, with their own stiffness:
`jointHertz` 60 and `jointDampingRatio` 2 by default, capped at half the substep rate. The
relax pass solves them rigid.

After restitution, two settle iterations move the bodies of every active joint directly:
the turn first, then any limit that is exceeded, then the point or the slide. Velocities
aren't touched.

## Consequences

The 12 link chain, total stretch at the end of 10 seconds, against Jolt with 10 velocity
and 2 position steps:

| Case | Twofold | Jolt |
| --- | --- | --- |
| Swinging, ball 209 times a link | 0.20 m | 2.22 m |
| Hanging still, 209 times | 0.17 m | 0.30 m |
| Hanging still, 21 times | 0.023 m | 0.046 m |

- A chain with a much heavier end still stretches some. A solver over the joint tree would
  fix that. Not worth it at these numbers.
- The settle pass knows nothing of contacts, so it can move a body slightly into
  something it touches. The next step's contacts push it back out.
- Scenes without joints step exactly as before.
