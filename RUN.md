# PReSto Temp12k run, 2026-10-08 pool

This repository runs the PReSto Temp12k reconstruction on a pinned LiPDverse data
pool. Everything that defines the run is committed here:

- `query_params.json` (`mode: bundle`): the pool bundle's URL and sha256,
  and its TSids. Pool, criteria and deduplication:
  [paleopresto/presto-recipes](https://github.com/paleopresto/presto-recipes/tree/pools/pools),
  release [pools-2026-10-08](https://github.com/paleopresto/presto-recipes/releases/tag/pools-2026-10-08).
- The algorithm's config, from presto-recipes `runs/presto-Temp12k/2026-10-08/`
  (or the template default where that directory has none).
- The workflow and container code, from the template's `pool-bundles` branch
  until that is merged upstream.

Committing `query_params.json` triggers the reconstruction workflow.
