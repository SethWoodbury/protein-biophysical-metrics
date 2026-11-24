# Protein biophysical metrics

Compute descriptive metrics for stored protein sequences and structures. Run
the same CLI on one sequence, a FASTA file, or a structure manifest.

```sh
pip install .
protein-biophysical-metrics --sequence MKT --output properties.csv
protein-biophysical-metrics --fasta proteins.fasta --groups composition,instability,optics --output properties.csv
protein-biophysical-metrics --manifest structures.csv --groups sap --workers 4 --output surface.csv
```

Manifest columns `id` and `sequence` are required. Optional `structure_path`
and `chain` select coordinates relative to the manifest and an author chain.
Structure metrics check exact sequence correspondence. CSV contains scalar
columns. Use `--jsonl PATH --details` to also get structured annotations and
their construct-relative coordinates. Plain `--jsonl PATH` writes scalar rows.

| Group | Measures |
|---|---|
| `composition` | Length, average molecular weight in Da, GRAVY, cysteine content, `methionine_count`, aromaticity, and aliphatic index |
| `charge` | Theoretical pI, signed and absolute net charge at `--ph`, with free unmodified termini |
| `instability` | Guruprasad sequence instability index |
| `optics` | Estimated 280 nm extinction coefficients for reduced and maximum disulfide-paired cases in M⁻¹ cm⁻¹ |
| `n-degron-stability-index` | Five-residue N-terminal windows and predicted Protein Stability Index (PSI) from N-FIVE. Lower PSI means a more destabilizing N-terminus. |
| `swi` | `swi_solubility_index`, the unitless sequence mean of the original SWI weights. Higher means higher model soluble propensity. |
| `site-flags` | `ng_count`, `dg_count`, `ns_count`, `polyproline_cluster_count`, `degron_tag_count`, `degron_motif_count`, and optional `protease_cut_site_count`. All are counts. |
| `aggrescan` | `aggrescan_na4vss`, normalized global profile score, and `aggrescan_thsar`, hotspot area per residue in model units. Higher means greater model aggregation propensity. |
| `helix` | `helix_end_charge_compatibility`, mean C-end charge minus mean N-end charge over eligible helices, in elementary-charge units per residue. |
| `sap` | Static SAP approximation, total, per-residue and side-chain SASA in Å², Cα radius of gyration in Å, and covered residues |
| `surface-chemistry` | `cys_sulfur_sasa_a2` (Cys SG) and `met_sulfur_sasa_a2` (Met SD) in Å², and unitless `hydrophobic_burial_fraction` in [0, 1]. |

The defaults are `composition,charge,instability`. Select other groups with
`--groups`. Pass `--context NAME:PREFIX:SUFFIX` repeatedly to evaluate explicit
constructs. For an untagged input, use `--context vector:MSG:SUFFIX`. For an
input already starting with the vector's `MSG`, use `--context vector:MSG+GGS:SUFFIX`
to insert `GGS` after it, or submit the full construct directly. Surface
metrics use the base row's coordinates. Supply a separate manifest row with
its own coordinates and sequence to evaluate a different proteoform.

Sequence descriptors characterize the stated construct and
the instability index is a sequence heuristic. The estimated pI is the zero
of the Bjellqvist sequence charge model. Use `--ph` between 0 and 14 for the
reported net charge. Cα radius of gyration uses the unweighted Cα centroid.
The aliphatic index describes composition, and estimated extinction depends
on the stated disulfide case.
[pI model](https://doi.org/10.1002/elps.11501401163) ·
[Aliphatic index](https://doi.org/10.1093/oxfordjournals.jbchem.a133168) ·
[Extinction coefficients](https://pmc.ncbi.nlm.nih.gov/articles/PMC2143013/)

### Sequence metrics and site details

```sh
protein-biophysical-metrics --fasta proteins.fasta --groups composition,swi,site-flags,aggrescan --output properties.csv
protein-biophysical-metrics --sequence QENLYFQSG --groups site-flags,helix \
  --proteases-used tev --intended-cuts 7 --protected-positions 1 \
  --mature-n-terminus --secondary-structure HHHHHHHHH \
  --output properties.csv --jsonl properties.jsonl --details
```

The chemical motifs `NG`, `DG`, and `NS` are counted separately. Details
record the sites, conditional `DP` cleavage liability, and possible `Q/E`
cyclization at a declared exposed mature N-end. Details also record classical
ClpS `F/L/W/Y` N-end residues and secondary `R/K` signals separately. Total
residue SASA in Å² is attached when already available from a selected
structure group. Overlapping `PPP` and `PPG` hits form one cluster. Details
keep the hits and longest Pro run. These are sequence and pausing annotations.
[Chemical liability](https://pubmed.ncbi.nlm.nih.gov/1649998/) ·
[Chemical stress](https://pubmed.ncbi.nlm.nih.gov/33475951/) ·
[Polyproline pausing](https://pubmed.ncbi.nlm.nih.gov/23239624/)

C-terminal flags apply to a free end of the complete construct, with ordinary
cytosolic BL21/BL21(DE3) as the biological scope. Complete native/LVA/AAV/ASV
ssrA tags take precedence over the longest matching shorter `LAA`, `VAA`,
`AA`, or `RRKKAI` motif. Details identify the match. Use `--blocked-c-terminus`
to declare a blocked input end. A context that supplies a new C-end uses the
free-end default. The match counts describe sequence context.
[ssrA variants](https://doi.org/10.1186/s13036-024-00420-9) ·
[C-terminal ClpX signals](https://pmc.ncbi.nlm.nih.gov/articles/PMC517522/) ·
[MuA recognition](https://doi.org/10.1016/j.jmb.2015.03.008)

`--proteases-used tev,hrv-3c,thrombin,enterokinase,factor-xa` accepts any
subset of the listed workflow enzymes. The optional count is blank/null
when none are selected. For a selected set it counts distinct unexpected
canonical cleavage bonds, including a zero for a complete scan with no such
matches. Details keep enzyme identities and intended matches. The registry
covers the stated recognition motifs, and off-target cleavage can occur
outside these patterns. Positions are 1-based. Intended cuts name the residue
immediately before the bond. Protected positions stay counted and are marked
in details. Contexts shift input positions for residues and bonds that remain.
Annotate a bond introduced by a prefix or suffix on its own full-construct row.
[TEV specificity](https://www.neb.com/en-us/products/p8112-tev-protease) ·
[Workflow enzyme specificities](https://b2b.sigmaaldrich.com/SE/en/technical-documents/technical-article/protein-biology/protein-purification/cleavage-enzymes-fusion-protein-purification)

For per-row declarations, optional manifest columns are `intended_cuts`,
`protected_positions`, `mature_n_terminus`, `free_c_terminus`, and
`secondary_structure`. Quote comma-separated position lists in CSV. Use
`true` or `false` for terminal declarations. Provide a processed sequence as
its own row to describe its actual ends. The helix group accepts a DSSP code
per residue and compares the first and last four residues of each contiguous
`H` run of at least eight residues. Asp/Glu contribute −1, Lys/Arg +1, and
all other residues 0. Higher values indicate greater compatibility under
this convention. Missing assignments or the absence of eligible helices
return null with a reason in details.

SWI uses an independently implemented arithmetic mean and the released
three-decimal residue weights credited to Bikash K. Bhandari, Chun Shen Lim,
and Paul P. Gardner (© 2019–2020). The weights are the third column of
`hydrop_flex_swi` in the authors' `functions.py` (SHA256
`bd1808337c372d2dcd2721fc28c2d25fd4c85b8e23ba52361d4135d7cfd9b48c`).
[SWI paper](https://doi.org/10.1093/bioinformatics/btaa578) ·
[Released weights](https://raw.githubusercontent.com/Gardner-BinfLab/TISIGNER-ReactJS/master/functions.py) ·
[TISIGNER software notice](https://github.com/Gardner-BinfLab/TISIGNER-ReactJS/blob/master/LICENSE.md)

AGGRESCAN follows the historical web-compatible psico conventions for
terminal padding, replicated edge windows, and the N-terminal hotspot area
correction. Pro breaks a hotspot. Na4vSS may be negative. Sequences of one
or two residues return null for both outputs. Details use 1-based inclusive
hotspot intervals. The pinned numerical adaptation and BSD notice ship in
`aggrescan.py`. [AGGRESCAN method](https://doi.org/10.1186/1471-2105-8-65) ·
[Pinned psico source](https://github.com/speleo3/pymol-psico/blob/f42cb1e4f565685b95fba63820b5f37db32b5954/psico/aggrescanning.py)

### N-degron-stability-index

From the repository root, initialize the pinned Kunjapur Lab submodule and
install the optional `nfive` dependencies. The Git LFS skip keeps the upstream
sequencing databases out of the checkout while still fetching the model.

```sh
GIT_LFS_SKIP_SMUDGE=1 git submodule update --init vendor/N-terminal-cluster-stability
python -m pip install -e '.[nfive]'
protein-biophysical-metrics --sequence ACDEFGH --groups N-degron-stability-index --context tagged:MSG: --output stability.csv
```

Each row reports `n_degron_sequence` and `n_degron_stability_index` for residues
1–5 of its construct. If the construct starts with `MSG`, it also reports
`n_degron_sequence_met_excised` and `n_degron_stability_index_met_excised`
for residues 2–6.
In the example, the base row scores `ACDEF` and the tagged row scores `MSGAC`
and `SGACD`. If the input already includes the tag, pass the construct directly.
Unavailable five-residue windows produce blank CSV cells and JSON `null`s.
The after-excision score describes the sequence that would be exposed if the
initiator Met is removed. The index is N-FIVE's predicted PSI in the authors'
*E. coli* reporter assay. Lower values indicate a more destabilizing N-terminus.
[N-FIVE source and model](https://github.com/KunjapurLab/N-terminal-cluster-stability) ·
[Sen et al.](https://doi.org/10.1038/s41467-026-79349-0)

For a non-editable install, set `--nfive-model` to the path of the trusted
`N-FIVE WT model.pkl` in the submodule. The Python API takes the same path as
`nfive_model` in `describe(sequence, groups=("n-degron-stability-index",), nfive_model=...)`.

N-FIVE batches unique windows and caches direct predictions in memory. To
generate the optional float32 `MSGXX`/`SGXXX` lookup from the trusted model
once, then reuse it across runs:

```sh
python -m protein_property_metrics.nfive --build-msg-lookup nfive-msg.npz
protein-biophysical-metrics --fasta tagged.fasta --groups n-degron-stability-index \
  --nfive-lookup nfive-msg.npz --output stability.csv
```

The lookup is checked against the alphabet and SHA256 of the selected model.
Other windows use direct prediction. In Python pass `nfive_lookup` to
`describe()` or use `prepare_nfive_scores()` for a batch of normalized constructs.

### Structure metrics

```sh
protein-biophysical-metrics --manifest structures.csv --groups sap,surface-chemistry,site-flags \
  --output surface.csv --jsonl surface.jsonl --details
```

The static SAP approximation uses atom-centered neighborhoods, accessible
side-chain surface and a glycine-centered hydrophobicity scale. The default
SAP radius is 5 Å. SASA uses Shrake–Rupley on the isolated canonical heavy-atom
chain with a 1.4 Å solvent probe and elemental van der Waals radii. The
default is 960 points per atom. Adjust sampling resolution with
`--sasa-points`. The positive residue total
depends on size. Its per-residue counterpart divides by the number of
covered residues. The signed maximum describes the highest scoring patch.
Compare values calculated using the same radius, SASA sampling and metric
definition. [SAP definition](https://doi.org/10.1073/pnas.0904191106) ·
[Reference side-chain surface areas](https://doi.org/10.1093/protein/15.8.659) ·
[Hydrophobicity scale](https://doi.org/10.1016/0003-2697(91)90045-U)

SAP and surface chemistry share a single atom-level SASA calculation. The
scope is the isolated canonical heavy-atom chain. Waters, metals, other
ligands and other chains are excluded from SASA. Modified residues in the
selected polymer, incomplete heavy-atom coverage and sequence mismatches
raise errors. Explicit PDB/mmCIF sulfur connections can appear in details
as chemical context. Sulfur totals are static exposure burdens. A covered
chain without a residue type gives zero for that type. Missing coordinates
give null. Hydrophobic burial uses glycine-centered positive Black–Mould
weights and clips exposed side-chain area divided by the reference area to
[0, 1]. Higher values mean more weighted hydrophobic burial. A chain with
zero total hydrophobic weight gives null.

The Python sequence API accepts `describe(sequence, groups=("swi", "site-flags"),
proteases_used=("tev",), intended_cuts=(7,), details=details)`, where `details`
is a dict populated in place. Select structure groups with
`describe_structure(path, expected_sequence=sequence, groups=("sap", "surface-chemistry"))`.

Run numerical and integration checks with `python -m unittest discover -s tests -v`.

### Measured CPU time

On `carnation` with one CPU thread and Python 3.12, 10,000 synthetic
150-residue sequences took 1.11–1.13 s for `composition,charge,instability,swi,site-flags,aggrescan`
with two workflow proteases enabled. All 8,400 distinct fixed-MSG N-FIVE
windows took 0.28–0.32 s to predict in one batch after model loading.
Repeating the batch took about 0.006 s from the in-process prediction cache
or 0.010 s from the loaded 35 KB lookup. These timings use in-memory sequences
and exclude output I/O. The first direct score took 3.1–4.5 s, including model
loading. First use of the lookup for all 8,400 windows took 0.049 s, including
file reads and model hash validation.
