# Production landscapes

EXP-066 writes the data for comparisons of global empirical and synthetic
landscape statistics here, grouped by RUN. This is the project path corresponding
to `processed.production_landscapes` in the scientist's request.

The queued job generates normal labproj outputs, verifies every SHA-256 in its
manifest, and atomically publishes the complete bundle to the durable cluster
counterpart:
`/tank/data/labproj/projects/project-adapt-env/data/processed/production_landscapes/<RUN-ID>`.
A failed or partial execution is retained in runtime outputs, not published as a
completed production dataset. Existing bundles are never silently overwritten.

The local counterpart is `data/processed/production_landscapes/<RUN-ID>`.
Use `scripts/collect_production_landscapes.py <RUN-ID>` to collect and verify it;
then `labproj data track data/processed/production_landscapes/<RUN-ID>` and
`labproj data push` to record and back up its DVC pointer. Runtime outputs and the
durable published bundle both retain the frozen config and exact source hashes.

Core products: `suitability.csv`, `selected_parameters.json`,
`confirmation_statistics.csv`, `figure_data.csv`, and per-assay raw cube arrays,
latent arrays, optimizer traces, mutation effects and observation sensitivities.
The figure table retains empirical panels and independent generator seeds.
Repeated observation draws and multiple cubes from one seed are not independent
landscapes. An empty high-order statistic for a two-site cube means unavailable.
