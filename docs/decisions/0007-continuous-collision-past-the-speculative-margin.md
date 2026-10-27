# 7. Continuous collision only past the speculative margin

Date: 2026-10-27

Status: Accepted

## Context

Speculative contacts stop a body moving at the velocity the step started with, so the
usual tunneling case is already handled. A body can still get through a thin shape when
its speed changes within the step: a hit from another body, or a very large force.

Options:

- Sweep every moving body every step. Expensive, and most sweeps find nothing.
- Mark some bodies as bullets and sweep only those. Requires the user to know in advance.
- More substeps. Costs everything and still fails at high enough speed.

## Decision

After the solve, a shape is fast when it moved further than its speculative contacts
covered plus half its inner extent. Only fast shapes are swept.

The sweep has two stages. First a sphere at the shape's core is swept with conservative
advancement. If the core never reaches the other shape, nothing tunneled and the body
keeps its pose. Only then is the full shape swept, and the body moves back to its first
touch. Dynamic bodies that hit each other move back together.

For a planar body, the inner extent ignores faces along the plane normal, so a thin layer
is judged by its width in the plane.

## Consequences

- No settings and no bullet flag. Nothing is swept in a normal step.
- Grazing contacts don't rewind bodies, so fast bodies don't stutter.
- A rewound body loses the rest of that step's motion.
