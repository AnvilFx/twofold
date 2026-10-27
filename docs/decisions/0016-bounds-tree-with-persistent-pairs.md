# 16. A bounds tree with persistent pairs

Date: 2026-10-27

Status: Accepted. Replaces [decision 4](0004-sort-and-sweep-rebuilt-every-step.md).

## Context

Sort and sweep sorted every shape every step, including static ones and sleeping piles.
That cost grows with the size of the world, not with how much of it moves. Cutouts spent
a quarter of its step there. Ray casts and the sweeps of continuous collision also tested
every shape.

## Decision

Every shape except half-spaces has a leaf in a dynamic AABB tree, balanced by rotations
and built by surface area cost as in Box2D. A leaf covers the shape's bounds grown by the
speculative distance and then by a 100 mm margin, and only moves in the tree when the
shape leaves it. Half-spaces have infinite bounds, so they stay in a short list of their
own and are tested against every moving shape.

The world keeps a sorted list of shape pairs whose leaves overlap. Only shapes whose leaf
was inserted or moved since the last step look for new pairs. Pairs whose leaves stopped
overlapping are dropped. Pairs between two bodies that are not dynamic aren't kept, so
changing a body's motion type sends its shapes looking again.

The candidates for the narrowphase are filtered from that list with exactly the old test:
bounds grown by the speculative distance and the body's travel, filters, and at least one
moving body. Ray casts, continuous collision and waking query the tree.

## Consequences

- Results are bit for bit what sort and sweep gave, on every scene and every ray.
- The broadphase costs about what moves. With 9,800 shapes and 40 awake bodies it went
  from 2.46 ms to 0.13 ms a step.
- The tree's shape depends on the order of inserts and removals, but nothing reads that
  order. Pairs are sorted by shape index, and rays try their hits in index order.
- More code than a sort, and a world holds a tree node per shape plus the pair list.
