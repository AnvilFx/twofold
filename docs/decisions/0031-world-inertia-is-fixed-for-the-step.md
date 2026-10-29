# 31. World inertia is fixed for the step

Date: 2026-10-29

Status: Accepted.

## Context

Each contact point direction keeps its lever arms, its velocity change per unit impulse and
its effective mass, built once per step. The world inverse inertia was rebuilt after every
position integration, and the angular responses with it, but the effective masses stayed as
they were at the start of the step. So a turned body took impulses sized for one inertia
and answered them through another.

For a round body that is harmless. A 1 cm by 1 cm by 1 m rod has about 5000 times more
inertia across it than along it, and a small turn changes its world inertia a lot. Dropped
flat from 1 m in 200 random orientations, 32 rods came back above their drop height, the
worst at 1.8e8 m with 2e10 rad/s of spin. CI found it when the macOS build, which fuses
multiply-adds in the tests, started a degenerate prism in a slightly different pose and it
was still in the air at the end of its drop test.

Options:

- Rebuild the effective masses along with the responses after every position integration.
  That costs the setup of every contact once per substep, and the lever arms would still be
  from the start of the step.
- Keep the world inverse inertia fixed for the whole step.

## Decision

Each body's world inverse inertia is refreshed in the step state and then kept until the
next step. Contact responses are no longer updated after each position integration, and
joint rows read the same inertia as the contacts. The settle pass after the substeps still
refreshes it for jointed bodies, since it moves them.

## Consequences

- None of the 200 rods comes back above its drop height, and none of 200 thin cylinders.
  Of 200 degenerate prisms, 1 still does, to 1.35 m, against 36 before.
- A 1 cm rod that lands flat at 3 m/s can still bounce up to 0.83 m with 4 substeps, from
  the two contact points across its width, which barely move it along but spin it hard about
  its length. 8 substeps settle it.
- A body that turns a lot within one step is solved with the inertia it had at the start.
- One pass over the active contacts per substep is gone.
