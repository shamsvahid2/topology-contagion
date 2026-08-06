# Networks

The empirical network files used by the paper are included here directly (carried over
from the repo's pre-existing `networks/` folder, which the scripts already reference via
a relative `NETWORK_DIR` — no path changes were needed):

- `unweighted/` / `weighted/G_*_conf.graphml` — RFID conference contact network (403
  nodes, 9,565 edges) — main dataset used alongside email in Figs. 2/3.
- `unweighted/` / `weighted/G_*_school.graphml` — Utah school network.
- `unweighted/` / `weighted/G_*_email_eu_modified.graphml` — EU email network (1,005
  nodes, 16,706 edges) — the headline dataset for Figs. 2, 3, 4.
- `weighted/G_weighted_hospital.graphml`, `weighted/G_weighted_work.graphml` — RFID
  hospital/workplace networks (SI robustness datasets).
- `unweighted/G_unweighted_barabasi_albert.graphml` — synthetic BA benchmark graph.

**Not included:** the large corpus of `G_*_smallworld_p*.graphml` Watts–Strogatz graphs
(hundreds of files, tens of MB) used for the SI small-world benchmark (SI Figs. S16/S17)
— these are regenerable via `notebooks/SI_small_world_generation.ipynb` /
`SI_small_world_generation_owen.ipynb` rather than checked in.

**Not included:** any simulation output or EPH-feature CSVs (the `result_*/` folders) —
per instruction, only code + the source network files are kept in this repo. Running the
scripts in `../contagion_simulation/` and `../eph_computation/` will regenerate these
locally into `result_*/` folders next to the scripts (create the folder first — the
scripts do not auto-create their output directory).
