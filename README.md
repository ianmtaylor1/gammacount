# gammacount
Distributions related to gamma renewal processes.

This R package contains distribution functions (density, probability, quantile, and random number generation) for gamma renewal processes. Distributions include counts for a standard and stationary gamma renewal process, the limiting excess life distribution, and subsequent renewal time distributions in the stationary process.

## Development

The `main` branch is for releases. The latest commit on the `main` branch should always be tagged with a release number. This is so that `devtools::install_github` should work with the defaults. Development takes place on other branches for feature updates, bugfixes, or documentation changes which are merged back into `main` when they are complete.

## Release numbering

This package's releases are numbered in a major, minor, bugfix system. The last number increments with bugfixes and documentation changes. Nothing about the interface to the code should change with an increase in this number. The second-to-last number is the minor release number. An increment of this number happens when features are added, but no existing features or interface aspects are changed or removed. Code written for an earlier minor release should still work on a later minor release. The first number is the major release number. This number will change when some existing code or feature has changed in a way that could make existing code incompatible with the new version.

Only the latest release is maintained at any one time. Use old versions at your own risk. No bugfixes will be released for old versions, so updating to fix a bug may require updating to a new minor or major release.

## Contributing

If you have a contribution, please send a pull request. If possible, make a separate branch based on the latest release tag to contain your changes. I am especially appreciative of:

* Improvements to the handling of edge cases (invalid parameter values, inputs not in the support of the distribution, etc.)
* Tests and increasing test coverage
* Documentation, especially clarifications
* Speed improvements

If you have a request for a change but cannot code it yourself, please let me know in a bug report or feature request.
