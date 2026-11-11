# Data license and attribution

ChEMBL database records are provided by EMBL-EBI under Creative Commons
Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0):
https://creativecommons.org/licenses/by-sa/3.0/

Source: https://www.ebi.ac.uk/chembl/ ; API:
https://www.ebi.ac.uk/chembl/api/data/

The exact database release, retrieval URLs/time and raw snapshot hashes are in
`raw/provenance.json`. Activity, molecule, assay and source-document identifiers
are retained. Curated tables are adaptations of this data and retain CC BY-SA 3.0.
The repository MIT license applies to original code, not these source records.

Changes: assay filtering, rejection of censored/flagged/multicomponent records,
unit conversion to pMIC, canonical identity grouping, replicate aggregation,
scaffold labels and training/test assignment. Original API payloads are archived
as gzip JSON. Consult individual source publications via their ChEMBL document
IDs for experimental provenance; no publication text or figures are redistributed.
