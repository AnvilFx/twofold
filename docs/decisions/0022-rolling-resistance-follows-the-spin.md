# 22. Rolling resistance follows the spin

Date: 2026-10-27

Status: Accepted

## Context

Exact spheres and cylinders roll forever on flat ground. Real ones stop because the
contact patch deforms. Engines usually model that as a torque against rolling, limited by
a coefficient times the normal force times the radius.

A torque solved per world axis would push a planar body about axes it is locked on, where
it has no effective mass.

## Decision

Rolling resistance is one angular impulse per contact, solved along the current relative
spin of the two bodies and then clamped to the limit as a vector. The limit is the mixed
coefficient times the total normal impulse times the mean height of the dynamic body's
center of mass above the contact points. The coefficient mixes by maximum.

## Consequences

- A ball rolling at v stops after 1.4 v^2 / (2 c g), a disc after 1.5 v^2 / (2 c g).
- A ball rests on a slope gentler than c and rolls down a steeper one.
- Spin about the contact normal is resisted too.
- With a coefficient of zero, the default, nothing changes and scenes step exactly as
  before.
