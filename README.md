# gammacount
The Gamma-count distribution, related distributions, and extensions

This R package contains functions for the gammacount (gc), first arrival time (ft), and gammacount with random start time (rst) distributions.

## Development

The `main` branch is for releases. The latest commit on the `main` branch should always be tagged with a release number. This is so that `devtools::install_github` should work with the defaults. Development takes place on other branches for feature updates, bugfixes, or documentation changes which are merged back into `main` when they are complete.

## Release numbering

This package's releases are numbered in a major, minor, bugfix system. The last number increments with bugfixes and documentation changes. Nothing about the interface to the code should change with an increase in this number. The second-to-last number is the minor release number. An increment of this number happens when features are added, but no existing features or interface aspects are changed or removed. Code written for an earlier minor release should still work on a later minor release. The first number is the major release number. This number will change when some existing code or feature has changed in a way that could make existing code incompatible with the new version.

Only the latest release is maintained at any one time. Use old versions at your own risk. No bugfixes will be released for old versions, so updating to fix a bug may require updating to a new minor or major release.
