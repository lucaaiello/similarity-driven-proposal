# Original records and current locations

`original_sources.zip` preserves the source/input/environment files that remained
from the completed application workflow before reorganization. Entry names are
the original recorded paths. The archive includes old wrappers and preparation
scripts for provenance; they are not the public run instructions. One
historical README containing personal login instructions is omitted; use
`../Cluster/README.md` instead. All retained archived entries preserve their
original bytes, and every source needed by recorded hash checks remains present.

`verify_sources.R` extracts only the named archived entries into a checked
temporary directory and computes their actual MD5 hashes. It verifies the
unchanged source receipts stored in the completed summary and input plans.
The saved records are not edited to manufacture agreement with renamed code.

`file_locations.csv` maps original paths to their current file or archive entry
and records their original SHA-256 values. `native_engine_receipt.rds` provides
the original native/reference checksums. The current native C++ and headers
match those original bytes. Historical packaging/environment files are retained
as records; use `../Cluster/` for the current directory layout and Slurm wrappers.

Location-map rows marked `metadata_redacted` identify saved metadata whose
personal paths were removed. Their original SHA-256 values refer to the
pre-redaction objects; current public input hashes are in
`../../reproduce/inputs_md5.csv`. Numerical values were not changed.

Live code for a new run is hashed independently by the runtime/audit gate. Any
new sampling still requires matching full checks in the target environment.
No claim of new sampler convergence follows from directory reorganization.
