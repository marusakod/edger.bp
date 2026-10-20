# Package website

The site uses the existing README, function documentation, and vignettes.
It is published at https://marusakod.github.io/edger.bp/.

The GitHub Actions workflow builds the site on pushes to `main`, validates
pull requests without deploying them, and can also be run manually. It uses
the audited R 4.5 dependency stack in `environment.yml` and installs the
current checkout before evaluating the examples.

For the initial deployment, a repository administrator must select **Settings >
Pages > Deploy from a branch**, then choose **gh-pages** and **/ (root)** after
the workflow has created that branch.

To build locally in the documented environment:

```bash
micromamba create -f pkgdown/environment.yml --strict-channel-priority
micromamba activate edger-bp-docs
R CMD INSTALL .
Rscript -e 'edger.bp::bp_check_versions(error = TRUE)'
Rscript -e 'pkgdown::build_site()'
```

On SHPC, run dependency installation and rendering through Slurm. Generated
HTML in `docs/` is ignored by Git and excluded from the source package.

To refresh the comparison table in the README, run the existing vignette with:

```bash
Rscript pkgdown/update-readme.R
```

This copies its six numerical comparisons and dependency versions into the
marked README section. Review and commit the resulting README change with
any changes to the vignette; do not edit its numerical results by hand.
