# Application sampler

- `dm_sampler.cpp`: native production transition engine, unchanged from the
  completed application. `dm_stable_numerics.h` and `joint_kernel.h` are its
  numerical/conditional-update headers, also unchanged.
- `runtime.R`: initializes persistent native contexts and one-thread workers,
  propagates local updates, and reconstructs contexts after serialization.
- `analysis_support.R`: fresh starts, the nine-edge exchange schedule, current
  source receipts and matching-runtime validation gate.
- `ladder_steps.R`: ladder iteration/checkpoint orchestration and safe cleanup
  of superseded checkpoints within a new run.

The common driver at `../run_analysis.R` selects a method/start task. Do not
invoke native helper exports as an alternative to the driver or bypass its
matching-source audit. New sampling requires full checks in the target runtime.

Original native-source generation receipts are in `../Provenance/`. Reference
engines are kept under `../Validation/Reference/` for comparison, rather than
presented as alternative application samplers. Moving and renaming the native
module did not modify its contents or numerical calculations.
