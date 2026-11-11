# Assessment of the completed fresh-start 50,000-iteration Slurm experiment

The main finding is a localized Prevotella exploration problem alongside much
stronger agreement of predictor PIPs outside that genus. The local sampler has
a measurable reproducibility advantage on most of the other genera, while
uniform flip has stronger cross-start agreement on the monitored Prevotella
fitted quantities at much lower observed computational cost.

The assessment uses all twelve downloaded summary records, each containing
25,000 retained cold draws. All source hashes matched the local audited source
graph, and all 21 full-window reported Rhat values were reproduced using
posterior::rhat (version 1.6.1) to absolute tolerance 1e-10. The downloaded summary
was unchanged after analysis. This assessment did not reread the remote full
checkpoints and did not advance a chain, compile, fit, or submit a job.

## Where agreement is strong

All PIP spreads of at least 0.5 across the four starts of any one method occur
in Prevotella. There are no such spreads among the other 29 genera.

For those other genera, the median of each genus's maximum cross-start PIP
spread is 0.03928 for similarity flip, 0.02100 for similarity local, and 0.09776
for uniform. The respective counts of genus-predictor entries with spreads of
at least 0.25 are 6, 3, and 9, out of 3,393 entries outside Prevotella.
Local has a lower genus-wise maximum spread than uniform in 28 of 29 genera;
Bacteroides is the exception. These are descriptive consistency measures,
not proof that local's estimates are closer to the unknown posterior PIPs.

Across methods, outside Prevotella only 13 of 3,393 entries have a difference
greater than 0.1 between method-averaged PIPs, and two exceed 0.25. Pooled
method means are exploratory where chains have not settled.

## The remaining subject4014 disagreement

The entries below are empirical means of fitted Prevotella fractions, in percent.

| Start | Similarity flip | Similarity local | Uniform |
|---|---:|---:|---:|
| Null | 0.00243 | 35.81 | 35.87 |
| 5% | 20.03 | 8.33 | 36.39 |
| 10% | 0.00305 | 0.00471 | 36.44 |
| 15% | 0.02182 | 36.03 | 43.20 |

Flip/null, flip/10%, flip/15%, and local/10% stay below the existing 5% fitted
threshold throughout all 25,000 retained iterations. This conclusion comes
from actual cold draws, independently of roundtrip counts.

The local/null ladder is an instructive counterexample to using zero full
roundtrips as evidence of failure: it has zero roundtrips but 5,746 end-of-sweep
5% crossings in the retained interval and a fitted mean close to uniform/null.
The local/15% ladder likewise has zero roundtrips, 2,054 threshold crossings,
and a fitted mean near 36%. The threshold describes a fitted-value region;
it is not an independently estimated posterior-mode boundary.

The local/5% ladder has 1,220 retained threshold crossings but substantial
changes in occupancy. Its fitted mean changes from 3.88% in the first retained
half to 12.78% in the second, then is only 0.54% over the final 5,000 draws.
Flip/5% changes from 14.26% to 25.80%, reaching 38.40% over the final 5,000.
Many crossings therefore do not alone settle the relative weight of regions.

Uniform samples both fitted regions in all four ladders and is the most
consistent on this subject, but the 15% ladder is higher and still drifting.
Full-window subject4014 Rhat is 1.305 for flip, 1.492 for local, and 1.036 for
uniform; diagnostic multichain bulk ESS is approximately 10, 7, and 109.
Uniform is a useful empirical comparator, not an established exact reference.

## What the retained windows add

Restricting to the last 12,500 draws does not remove the informed subject4014
gap: Rhat is approximately 1.538 for flip and 1.532 for local. Several local
diagnostics improve in the later window, including model size (1.060 to1.039)
and subject3026 (1.079 to1.048), but Prevotella model size does not improve
(1.317 to1.444). Uniform overall model-size Rhat also increases from1.025
to1.094 in that window. There is no justification here for retrospectively
selecting a shorter retained window as a cure for the disagreements.

Subject3026 shows a more favorable late trajectory for the informed methods.
For example, flip/15% changes from a fitted mean of5.48% in the first retained
half to31.02% in the second and40.19% in the final5,000 draws. This is evidence
of continued movement, although it is also evidence that the original retained
window includes a material transient.

## Associations with consistent empirical support

Seven genus-predictor entries have PIP at least0.5 in every one of the12 chains:
Bacteroides/Maltose, Bacteroides/Taurine, Bacteroides/Catechin,
Barnesiella/Calories, Prevotella/Arachidonic fatty acid,
Phascolarctobacterium/Iodine, and Parasutterella/Alcohol.
See all12_pip_ge_half.csv for exact chain minima and method means. Only two
entries have PIP at least0.8 in every chain, and none reach0.9 in every chain.
These thresholds describe cross-chain empirical support, not causal evidence
or a substitute for convergence assessment.

## Cost and presentation

Median wall time per ladder, including startup, is50.03h for similarity flip,
53.26h for similarity local, and3.44h for uniform. Total allocated core-hours
are approximately1,996,2,154,and137, respectively. These are descriptive
costs; load and visited states differ. The present experiment does not establish
an overall efficiency advantage for the informed methods on this dataset.

The real-data study can be presented as a transparent stress test, reporting
the strong agreement across most genera, the local sampler's reproducibility
advantages outside Prevotella, and the unresolved Prevotella differences.
An across-the-board convergence or superiority claim is not supported.

The next useful assessment is of the interpretable associations and correlated
predictor groups using saved results, with separate reporting of fitted-value
agreement and individual-variable selection. Further large sampler runs need
a specific hypothesis grounded in these results, rather than iteration count
alone.

## Scope of the downloaded data

The summary contains retained scalar/fitted traces and all-window PIPs, but
not every retained indicator matrix. Consequently this assessment cannot
measure first-versus-second-half PIP drift or indicator-specific ESS. Those
would require additional read-only postprocessing of the saved checkpoints.
The raw data and all scientific sampler sources/results were preserved.

Detailed tables are in report.txt and the accompanying CSV files. The five-page
retained_trace_assessment.pdf contains retained fitted/scalar traces, fitted
histograms and quarter-window comparisons. Faint trace lines are raw draws;
strong lines are125-draw block means.
