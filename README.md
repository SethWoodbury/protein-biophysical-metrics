# Protein biophysical metrics

Compute descriptive metrics for stored protein sequences and structures. Run
the same CLI on one sequence, a FASTA file, or a structure manifest.

```sh
pip install .
protein-biophysical-metrics --sequence MKT --output properties.csv
protein-biophysical-metrics --fasta proteins.fasta --groups composition,instability,optics --output properties.csv
protein-biophysical-metrics --manifest structures.csv --groups sap --workers 4 --output surface.csv
```

Manifest columns are `id,sequence,structure_path,chain`. The last two are
optional. Structure paths are relative to the manifest. Select an author
chain if a structure has multiple protein chains. The command checks exact
sequence correspondence when a structure is supplied for SAP. CSV is the
primary output. Use `--jsonl PATH` for a second copy as JSON Lines.

| Group | Measures |
|---|---|
| `composition` | Length, average molecular weight in Da, GRAVY, cysteine content, aromaticity, and aliphatic index |
| `charge` | Theoretical pI, signed and absolute net charge at `--ph`, with free unmodified termini |
| `instability` | Guruprasad sequence instability index |
| `optics` | Estimated 280 nm extinction coefficients for reduced and maximum disulfide-paired cases in M⁻¹ cm⁻¹ |
| `n-degron-stability-index` | Five-residue N-terminal windows and predicted Protein Stability Index (PSI) from N-FIVE. Lower PSI means a more destabilizing N-terminus. |
| `sap` | Static SAP approximation, total, per-residue and side-chain SASA in Å², Cα radius of gyration in Å, and covered residues |

Choose metric groups with `--groups`. Pass `--context NAME:PREFIX:SUFFIX`
to calculate sequence properties with explicit added residues. This flag is
repeatable. The base row uses the input sequence. Each context row uses the
concatenated sequence. Surface properties describe the supplied structure
on the base row. Sequence descriptors characterize the stated construct and
the instability index is a sequence heuristic. The estimated pI is the zero
of the Bjellqvist sequence charge model. Use `--ph` between 0 and 14 for the
reported net charge. Cα radius of gyration uses the unweighted Cα centroid.
The aliphatic index describes composition, and estimated extinction depends
on the stated disulfide case.
[pI model](https://doi.org/10.1002/elps.11501401163) ·
[Aliphatic index](https://doi.org/10.1093/oxfordjournals.jbchem.a133168) ·
[Extinction coefficients](https://pmc.ncbi.nlm.nih.gov/articles/PMC2143013/)

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

### Static SAP

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

Run the small numerical checks with `python -m unittest discover -s tests -v`
after installing the package. They cover theoretical pI and charge, optical
estimates, N-FIVE windows and encoding, analytical sphere areas, and PDB/mmCIF parity.
