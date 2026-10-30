# 33. SAT tests run four candidates at a time

Date: 2026-10-30

Status: Accepted.

## Context

With the contact solve batched, collision took 35% to 65% of a step in Rain, Hulls, Sweeper
and Meshes. Empty pairs cost little: 0.03 to 0.06 ms per frame in each of those scenes. The
cost is in pairs that touch. A touching hull pair took about 3 µs, a mesh pair about 6.5 µs.
Fast pairs are collided again within the step, which roughly doubles that.

In Hulls the time went to the edge query (21% of a step), the face queries (14%), the GJK
check that follows a pruned edge query (11%), finding the face toward a direction (7%) and
the box test that picks the edges near the overlap (6%). The face and vertex searches read
the best candidate's dot product again on every step of the loop.

Options:

- Warm start GJK from the line between the centers, or from the first search's normal for
  the core search.
- Stop GJK at an absolute gap tolerance instead of a relative one.
- Cache a separating axis for each pair that found no contact.
- Make the tests that decide the result cheaper without changing what they decide.

## Decision

The tests stay exact.

- The Gauss map tests of an edge query take four edges of B at a time with SSE2 or NEON,
  and the pairs that pass are finished one at a time in the same order as before. A
  triangle's edge query does the same.
- A face query scans the vertices of B four at a time. The first block with a vertex at or
  below the bound ends the scan, and the first such vertex in order is the result, as in
  the scan one vertex at a time. Otherwise the least of every lane is kept with `minf`
  semantics, which don't depend on order.
- The search for the face toward a direction and for the deepest vertex keep the best dot
  product instead of recomputing it.
- The edge box test marks each vertex with the sides of the box it lies past, and drops an
  edge when its two ends share a side.
- Where only the question "further apart than a limit?" matters, GJK stops once its upper
  bound is within the limit or its lower bound is past it. Its answer is the one the full
  search would give.

Lanes take the same float operations as the scalar test, so `TWOFOLD_SIMD=OFF` gives the
same bits, and every scene keeps its hash.

## Consequences

- Whole steps on one core, 12 seeds of 600 frames, best of three, before and after:
  - With sleep off: Hulls 0.529 to 0.470 ms, Rain 0.687 to 0.653 ms, Sweeper 0.328 to
    0.308 ms, Meshes 1.414 to 1.384 ms. Other scenes are within noise.
  - With sleep on: Hulls 0.437 to 0.326 ms, Meshes 1.388 to 1.319 ms, Sweeper 0.291 to
    0.275 ms.
- Rejected:
  - The GJK warm starts barely changed the iteration count (11.0 to 10.1 in Rain). Its
    long searches come from curved shapes near touching, not from a poor start.
  - An absolute GJK tolerance of 0.1 mm saved 0.05 ms in Rain and 0.07 ms in Meshes, but
    it lets the normal drift where the gap is flat near its peak. The 0.5 mm certificate
    was dropped for the same reason in an earlier change.
  - A separating axis cache could save at most the 0.03 to 0.06 ms that empty pairs cost.
- A pair with more than a few edges near the overlap still pays for every edge pair it
  tests. The second collision of fast pairs within the step is now the larger cost.
