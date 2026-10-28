# 27. Swing-twist rows follow the exact rates

Date: 2026-10-28

Status: Accepted

## Context

A ragdoll needs a shoulder or hip: a point joint whose limb stays inside a cone, often
wider one way than the other, and whose twist about the limb stays within a range. Doors
and drawers also need friction so they stay where they are left, and a game needs to know
how hard a joint pulls and to let it break.

Swing and twist come from splitting body B's turn in A's joint frame into a twist about
the frame's x axis and a swing that tilts it, as Jolt does. The cone is an ellipse in the
sines of the half angles, so it stays smooth up to a half turn.

The first limit rows turned about B's axis for twist and about the swing axis for the
cone. Those are where the angles change only while the limb hangs straight. Held at the
twist limit by a steady torque, a swung limb tilted out of its cone, since each limit row
pushed the other's angle.

Options:

- Rows about those simple axes, with more iterations. Slower, and the tilt stays.
- Swing as two separate angle limits, a box instead of an ellipse. Simpler rates, but
  corners that catch the limb.
- Rows along the exact rates: the relative spin that changes the twist, and the one that
  changes the room left in the cone.

## Decision

Twist is limited by two one-sided rows along `(xA + xB) / (1 + xA . xB)`, where xA and xB
are the joint axis on each body. That is the relative spin that changes the twist and
nothing else. The row is skipped when the swing is close to a half turn, where twist has
no meaning.

The cone is one one-sided row whose error is the angle left before the ellipse. Its rate
is the swing direction plus a term along `xB x k`, scaled by how fast the ellipse's reach
changes as the swing direction turns about the axis. k is the swing direction.

Both are anchored turns, like the hinge's limits in
[decision 18](0018-limits-and-motors-move-the-whole-joint.md). Rows go friction, twist,
swing, point. The settle pass moves any exceeded twist, then any exceeded swing, then the
anchors.

Friction on a hinge or slider is a motor row with speed zero, solved before the motor and
clamped to the friction torque or force times the substep. On a swing-twist joint it is a
3D turn row clamped by its length. Friction only adds to the warm start when it is set,
so a joint without it steps exactly as before.

The reaction is the impulse each row ends a substep's relax pass with, summed over the
step and divided by it: a force on body B and a torque about anchor B. The last substep
alone read about 2% high for a hanging body.

A joint breaks after the substeps and before the settle pass, when its reaction is over
`breakForce` or `breakTorque`. Its impulses are dropped and it leaves the active list.
`brokenJoints` lists the joints that broke in the last step, and the handle stays valid
until it is destroyed.

## Consequences

- A limb held at its twist limit under 0.5 N m stays within its cone. Under 2 N m on a
  light limb it slowly winds about the line between the two axes, which is free with the
  twist at its limit.
- A falling ragdoll stays together from 1 m and 10 m and in a pile of six. At rest the
  limits hold within 0.002 rad. In a frictionless pile, contacts can push a limb
  0.17 rad past its cone and 0.33 rad past its twist for a moment.
- A planar chain of swing-twist joints moves like the same chain of limited hinges,
  within 0.02 mm.
- Breaking is decided once a step from the summed reaction, so a joint breaks in the same
  step on every run and in every world copy.
- A broken joint keeps its collision filter. Bodies it joined still don't collide unless
  `collideConnected` is set or the joint is destroyed.
- Scenes without the new features keep their hashes.
