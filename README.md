# Protein property metrics

Compute descriptive metrics for stored protein sequences and structures. Run
the same CLI on one sequence, a FASTA file, or a structure manifest.

```sh
pip install .
protein-property-metrics --sequence MKT --output properties.csv
protein-property-metrics --fasta proteins.fasta --groups composition,instability --output properties.csv
protein-property-metrics --manifest structures.csv --groups sap --workers 4 --output surface.csv
```

Manifest columns are `id,sequence,structure_path,chain`. The last two are
optional. Structure paths are relative to the manifest. Select an author
chain if a structure has multiple protein chains. The command checks exact
sequence correspondence when a structure is supplied for SAP. CSV is the
primary output. Use `--jsonl PATH` for a second copy as JSON Lines.

| Group | Measures |
|---|---|
| `composition` | Length, average molecular weight in Da, GRAVY, cysteine content, aromaticity |
| `charge` | Estimated pI, signed and absolute net charge at `--ph`, with free termini |
| `instability` | Guruprasad sequence instability index |
| `sap` | Static SAP approximation, solvent-accessible surface area, covered residues |

Choose metric groups with `--groups`. Pass `--context NAME:PREFIX:SUFFIX`
to calculate sequence properties with explicit added residues. This flag is
repeatable. The base row uses the input sequence. Each context row uses the
concatenated sequence. Surface properties describe the supplied structure
on the base row. Sequence descriptors characterize the stated construct and
the instability index is a sequence heuristic.

The static SAP approximation uses atom-centered neighborhoods, accessible
side-chain surface and a glycine-centered hydrophobicity scale. The default
radius is 5 Å with 100 SASA points per atom. The positive residue total
depends on size. Its per-residue counterpart divides by the number of
covered residues. The signed maximum describes the highest scoring patch.
Compare values calculated using the same radius, SASA sampling and metric
definition. [SAP definition](https://doi.org/10.1073/pnas.0904191106) ·
[Reference side-chain surface areas](https://doi.org/10.1093/protein/15.8.659) ·
[Hydrophobicity scale](https://doi.org/10.1016/0003-2697(91)90045-U)
