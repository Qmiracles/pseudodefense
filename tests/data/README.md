# Test fixtures

Real outputs reduced for tests (no genomes are stored in the repo).

- `defensefinder/GCF_012971705.1/`: DefenseFinder 3.0.0 (models 3.1.0) systems and
  genes tables for RefSeq assembly GCF_012971705.1 (*P. aeruginosa* LYT4), run on
  pyrodigal 3.7.1 proteins with `--db-type gembase`. Unmodified; note that the
  systems table contains one exact duplicate line produced by DefenseFinder.
- `genes/GCF_012971705.1.gff`: pyrodigal GFF3 lines for the proteins in those systems only.
