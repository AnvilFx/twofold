# 18. Limits and motors move the whole joint

Date: 2026-10-27

Status: Accepted

## Context

A door swinging into its hinge limit bounced back with 14% of the speed it hit with. A
slider whose anchor sits away from the body's center bounced 38% or more.

The limit row only knew how the body turns about its own center. The joint makes it turn
about the anchor, so the point row then had to drag the center around. With Gauss-Seidel
the two rows trade impulses, and the warm started point impulse came back out as spin in
the wrong direction.

Tried and ruled out:

- No warm starting for the limits. No change.
- No settle pass. No change.
- Limits after the point row, the order Jolt uses. Worse.
- Softer joints. Less bounce, but chains stretch about four times as far.

## Decision

The hinge's motor and limit rows apply an angular impulse about the axis together with the
point impulse per unit of turn that keeps the anchors together. Their effective mass is
that of turning about the anchor.

The slider's slide, motor and limit rows include the turn impulse that keeps the relative
rotation fixed.

Rows go in the order Box2D uses: motor, limits, turn, point.

## Consequences

- The door and the offset slider stop at their limits within 2 mm or 0.002 rad, and are
  still a second later. A test covers each.
- Each of these rows costs a 3x3 pseudo-inverse per joint per pass.
- The composite row is exact for the joint on its own. A body that is also in contacts or
  other joints is still coupled to them by Gauss-Seidel.
- Hinges and sliders without limits or motors solve as before.
