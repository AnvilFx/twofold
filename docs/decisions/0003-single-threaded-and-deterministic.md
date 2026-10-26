# 3. Single threaded and deterministic

Date: 2026-10-26

Status: Accepted

## Context

Users scrub back and forth through a simulation and expect the same result each time.
Debugging is also far easier when a run can be replayed exactly.

## Decision

The step runs on one thread. Nothing depends on memory addresses or hash order: candidate
pairs and contacts are sorted by shape index, and bodies are visited by index. `World` keeps
all its state by value, so copying it gives a snapshot that steps exactly like the original.

## Consequences

- The same calls in the same order give bit-identical results with the same build.
- Rewinding is a copy, and replay can be checked frame by frame. The viewer does this.
- No multithreading. If it comes later, it has to keep the same results, for example by
  solving islands in parallel and merging them in a fixed order.
- Different compilers, flags or platforms can still differ, since floating point results
  depend on them.
