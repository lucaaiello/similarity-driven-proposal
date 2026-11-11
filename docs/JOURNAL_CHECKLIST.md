# Reproducibility release checklist

The package supplies executable code, data, saved results, a software/environment
record and a figure/table map. This follows JCGS's emphasis on enabling reviewers
to check reported results; see the
[ASA interview with its editor](https://magazine.amstat.org/blog/2017/10/01/dianne-cook/).

The publisher's current instructions page returned HTTP 403 to the automated
browser during preparation. Check the live
[JCGS instructions](https://www.tandfonline.com/action/authorSubmission?show=instructions&journalCode=ucgs20)
at upload. These checks do not constitute journal certification.

## Implemented locally

- One saved-results command for all figures and computational table summaries.
- Numerical checks against all 16 computational tables; three definition tables
  and the model diagram explicitly identified as static assets.
- Input checksums, software versions, timings, logs and completion status.
- Optional recomputation of parameter diagnostics from full retained fits.
- Data provenance, computational settings and optional sampling instructions.
- Verified journal archive excluding backups, Git history and correspondence.
- Explicit Git LFS rule for the large application result; journal archive uses
  actual objects rather than LFS pointers.

## Before public release/submission

- Select the code license with coauthors; intentionally undecided for now.
- Confirm terms for redistributing the supplied third-party data.
- Publish the prepared revision and check a fresh download.
- Archive an immutable release, preferably with a persistent DOI; insert the
  actual identifier in the manuscript availability statement when assigned.
- Upload the complete code/data archive as supplementary material, following
  the portal's actual size and anonymity rules. A GitHub link complements it.
- Ensure availability statements match the final location/access conditions.

Funding and other submission declarations remain separate author responsibilities.
The scripts do not publish or upload anything automatically.
