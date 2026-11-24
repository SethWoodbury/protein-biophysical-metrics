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
estimates, analytical sphere areas, and PDB/mmCIF parity.
