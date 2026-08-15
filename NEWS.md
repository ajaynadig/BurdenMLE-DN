# BurdenMLEDN 0.1.1

- Added a prepared gene posterior sampler that draws continuous effects from
  the Poisson likelihood within each selected fitted mixture component.
- Standardized full, bootstrap, and null fit status across MixSQP and EM;
  uncertainty now distinguishes usable nonconvergence from unusable fits.
- Made bootstrap samples gene-identified and added independent, replayable
  bootstrap and null seeds without changing caller RNG state.
- Made bootstrap reconstruction and confidence-interval shapes robust for
  one or many annotation strata and replicates.
- Disabled optional null simulations by default and enabled gene posterior
  summaries by default, matching routine substantive use.
- Added helpers for prevalence-derived and conservatively wide effect grids.
- Added count-free LOEUF strata and mutation-rate reference resources.
- Expanded documentation of output fields and manuscript-specific settings.

# BurdenMLEDN 0.1.0

- Initial public package structure.
- MixSQP fitting, gene-level bootstrap uncertainty, mutational variance,
  effective penetrance, posterior summaries, and null simulation.
- Synthetic worked example and manuscript-reproduction manifest.
