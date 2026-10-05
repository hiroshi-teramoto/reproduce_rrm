# GRRM preprocessing corrections

Initial validation: 2026-09-25. D*h review correction and full-batch audit: 2026-10-02.
Originally developed on fork `main` at `c46c62f`; the first PR revision was
`31c15f8`. The final changes from
[fork PR #25](https://github.com/mochimochi-brian/reproduce_rrm/pull/25)
(`d517c5e`) are ported onto upstream `main` at `d64568d`.

## Upstream port validation (2026-10-05)

Both regression commands below passed on qc4 with the `rrm-recon` Python
environment. The preprocessor, extractor, and both tests are byte-identical
to fork PR #25 at `d517c5e`; Python syntax checks and `git diff --check` passed.
The README additions were inserted into the upstream README without its
fork-specific surrounding text.

The calculation results below are retained evidence from the original fork
validation. Those CNPI runs used the fork's fast GAP implementation; this
port has no dependency on that implementation or the separate speedup PR.
No production CNPI calculations were rerun for the port.

## Behavior

`rrm_reconstruction_v18.py` treats linear `C*v` structures as achiral for EQ
and TS. In pymatgen 2024.8.9, the finite operation lists are only `[E]` for
`C*v` and `[E, i]` for `D*h`. Missing mirror operations caused the original
determinant-only chirality check to add spurious `C*v` inversion copies.

The same incomplete lists also omitted proper atom permutations for `D*h`:
filtering out inversion discarded the endpoint exchange realizable by a
180-degree rotation perpendicular to the molecular axis. The shared
`rotation_permutations()` now includes all listed operations for `C*v` and
`D*h` when constructing both `ur` (EQ) and `urt` (TS). A mirror plane through
a linear molecule fixes every atom, so composing it with an improper operation
produces the same permutation by a proper rotation. Nonlinear structures
continue to use only operations with positive determinant.

`extract_largest_component.py` explicitly selects the unique largest original
EQ/TS component by EQ count. It renumbers EQ and TS independently, updates
connections, preserves loops and parallel edges, copies retained TS logs byte
for byte, and records old IDs and excluded components in `mapping.json`.
Ties, invalid connections or IDs, missing retained TS logs, and existing output
paths cause errors. Failed writes clean up the newly reserved output directory;
parent directories may remain and are not removed because they may be shared.

See the [README](../../README.md#extracting-the-largest-connected-component)
for extraction and preprocessing commands.

## Regression checks and scope audit

PBS `756495.qc-cluster3.local` ran both tests and audited all EQ/TS geometries
in the original 29-input batch, using one CPU on `compute-52-0.local`.
The job moved from the full `lightq` to `welterq` before execution.

```sh
python3 tests/test_chirality.py
python3 tests/test_extract_largest_component.py
```

- Chirality and rotation permutations: HCN/CCO (`C*v`), OCO (`D*h`), bent
  water (`C2v`), and chiral CHFClBr (`C1`), including their inversions.
- OCO must have exactly identity and the `(1,3)` endpoint exchange.
- Nonlinear methane (`Td`) must have exactly 12 even permutations of its
  four H atoms; its 12 improper permutations must remain excluded.
- Extraction: renumbering, mapping, byte preservation, loops, parallel edges,
  invalid input, ties, missing logs, output protection, EOF handling, and
  connected/single-EQ inputs all passed.

The audit analyzed coordinates with the same pymatgen settings as preprocessing;
it did not select inputs only by saved `C*v` warnings. Six inputs contain
`D*h` and change under this correction: CO2 twice, CS2 twice, C2H2 EQ0, and
CH2 TS0. The other 23 inputs have unchanged rotation permutation sets.

## Corrected CNPI outputs

The first 14 reruns (`652126`–`652139`, 2026-09-25) fixed `C*v` inversion
copies, but their CO2/CS2 counts of **5 vertices / 4 edges were still wrong**.
The six `D*h` inputs were rerun in PBS `756496`–`756501` on 2026-10-02,
using `welterq`, fast, `GAP_WORKERS=1`, one CPU and 16 GB per job. All six passed.
Previous outputs are retained. The table gives the latest validated output;
unchanged `C*v`-only inputs retain their September 25 result.

| Input under the calculation directory | Latest PBS job | Before linear-symmetry fixes: vertices / edges | Latest vertices / edges |
|---|---:|---:|---:|
| `Catalyst/AgCO/AgCO_ADDF` | 652126 | Preprocessing disconnected | 2 / 1 |
| `Catalyst/AgH2/AgH2_ADDF` | 652127 | 1 / 4 | 1 / 2 |
| `Catalyst/AgOH/AgOH_ADDF` | 652128 | Preprocessing disconnected | 1 / 1 |
| `Catalyst/AuCO/AuCO_ADDF` | 652129 | Preprocessing disconnected | 2 / 1 |
| `Organic/C2H2/ADDF/wB97XD_6-31G*/C2H2_ADDF` | 756496 | Vertex lookup failed | 4 / 4 |
| `Organic/C2O/AFIR/B3LYP_6-31Gdp/C2O_AFIR` | 652130 | 5 / 4 | 3 / 2 |
| `Organic/CH2/ADDF/B3LYP_6-31Gdp/CH2_ADDF` | 756497 | 1 / 2 | 1 / 1 |
| `Organic/CNH/ADDF/wB97XD_6-31G*/CNH_ADDF` | 652131 | Preprocessing disconnected | 2 / 1 |
| `Organic/CO2/ADDF/B3LYP_6-31Gdp/CO2_ADDF` | 756498 | 7 / 6 | 4 / 4 |
| `Organic/CO2/ADDF/PBEPBE_6-31Gdp/CO2_ADDF` | 756499 | 7 / 6 | 4 / 4 |
| `Organic/COS/ADDF/B3LYP_6-31Gdp/COS_ADDF` | 652134 | 6 / 5 | 4 / 3 |
| `Organic/CS2/ADDF/B3LYP_6-31G2dp/CS2_ADDF` | 756500 | 7 / 6 | 4 / 4 |
| `Organic/CS2/ADDF/B3LYP_6-31Gdp/CS2_ADDF` | 756501 | 7 / 6 | 4 / 4 |
| `Organic/OCS/ADDF/B3LYP_6-31G2dp/OCS_ADDF` | 652137 | 7 / 6 | 4 / 3 |
| `Organic/OCS/ADDF/B3LYP_6-31Gdp/OCS_ADDF` | 652138 | 6 / 5 | 4 / 3 |
| `Organic/OCS/ADDF/UB3LYP_6-31Gdp/OCS_ADDF` | 652139 | 4 / 5 | 3 / 3 |

The new runs check that each affected `ur`/`urt` has order 2 and contains the
endpoint-exchange permutation before expanding it. They also passed checks
for unwanted inversion copies, preprocessing connectivity, vertex mapping,
group-index counts, output hashes and degree consistency. No Pechukas violation
was reported. C2H2 now completes instead of failing vertex lookup; CH2 loses
one spurious edge (2 -> 1) after restoring its TS rotational symmetry.

Group-index counts and degree consistency alone are insufficient to establish
correct symmetry: they use the same potentially wrong input groups. The new
known-geometry permutation tests supply the missing independent expectations.

Where needed, the reruns retained the existing input-level normalization
`Group([])` -> `Group([()])`, preserving the unnormalized input and both hashes.
That normalization is not a preprocessor source change in this PR.

## Pt5 extraction validation

PBS `651721` on `compute-46-0.local` (2026-09-21, one CPU) verified:

- 43 EQ / 76 TS / 6 components -> 38 EQ / 75 TS / 1 component.
- Excluded isolated EQs 0, 12, 29, 41 and 42, including TS7 at EQ12.
- Correct numbering/connections; unchanged remaining list bytes and all 75
  retained individual TS logs; unchanged hashes of all source logs.

The extractor logic is unchanged since that run; the review removed a stale
docstring reference and clarified parent-directory cleanup. Production CNPI
expansion of the six disconnected metal inputs remains pending.

## Provenance and raw artifacts

| Current source | SHA-256 |
|---|---|
| `rrm_reconstruction_v18.py` | `b3862104bd09aa1b869bb6900af9b700aaf4b21ca2ac77484e997802c4f1da13` |
| `extract_largest_component.py` | `b6eda896b2c8fec10c8ec46d5e608dd97a2d81d09591d93d2df9b4f2020fbf88` |

The current preprocessor hash matches the six October 2 calculation snapshots.
The September 25 runs used preprocessor hash
`0a2b9d6b667aa2b219a9e9aca3d93f6a97be399033ee4ab8833ea931c18e7c83`;
Pt5 extraction used `34c06d27549eeda65a276b8bee263c04dd2e053291510b39200a6573217b8a8a`.

Raw files remain outside Git, under the corresponding calculation directories:

- D*h correction: `CNPI/retry_linear_symmetry_20261002/` (snapshots, patches,
  logs, group-membership checks, and the `current/` output bundle).
- C*v-only results: `CNPI/retry_chirality_20260925/`.
- Full audit/tests: `Catalyst/AgCO/CNPI/pr_linear_review_20261002/`.
- Pt5 extraction: `Metal/Pt5/CNPI/extractor_validation_20260921/`.

Completed PBS history is disabled; success is recorded in compute-node status
files and validation logs. Source snapshots and input hashes are retained.
