## v0.3.1 - 2024-11-15

<!-- Release notes generated using configuration in .github/release.yml at main -->
### What's Changed

#### Bug Fixes

* Fix compatibility of interactive viewer with recent versions of Matplotlib by @astrofrog in https://github.com/dendrograms/astrodendro/pull/202

#### Other Changes

* Bump actions/checkout from 4.2.0 to 4.2.2 in /.github/workflows in the actions group by @dependabot in https://github.com/dendrograms/astrodendro/pull/201

### New Contributors

* @dependabot made their first contribution in https://github.com/dendrograms/astrodendro/pull/201

**Full Changelog**: https://github.com/dendrograms/astrodendro/compare/v0.3.0...v0.3.1

## v0.3.0 - 2024-11-15

<!-- Release notes generated using configuration in .github/release.yml at main -->
### What's Changed

* update "Iterable" import for Python versions > 3.9 by @nbrunett in https://github.com/dendrograms/astrodendro/pull/184
* Fix viewer in Python 3 by @indebetouw in https://github.com/dendrograms/astrodendro/pull/181
* Remove np.int by @keflavich in https://github.com/dendrograms/astrodendro/pull/179
* h5py removed .value after v3 by @indebetouw in https://github.com/dendrograms/astrodendro/pull/182
* Update package infrastructure by @astrofrog in https://github.com/dendrograms/astrodendro/pull/186
* Update plot.py to pass subtree to get_lines by @tonywong94 in https://github.com/dendrograms/astrodendro/pull/168
* Fix for changes between astropy 5.2.2 and 5.3.1 by @ajrigby in https://github.com/dendrograms/astrodendro/pull/192
* Fix tests/continuous integration by @astrofrog in https://github.com/dendrograms/astrodendro/pull/193
* WCSAxes import fix by @Parkerwise in https://github.com/dendrograms/astrodendro/pull/198
* Update infrastructure and fix compatibility with Numpy 2.0 by @astrofrog in https://github.com/dendrograms/astrodendro/pull/200

### New Contributors

* @nbrunett made their first contribution in https://github.com/dendrograms/astrodendro/pull/184
* @indebetouw made their first contribution in https://github.com/dendrograms/astrodendro/pull/181
* @tonywong94 made their first contribution in https://github.com/dendrograms/astrodendro/pull/168
* @ajrigby made their first contribution in https://github.com/dendrograms/astrodendro/pull/192
* @Parkerwise made their first contribution in https://github.com/dendrograms/astrodendro/pull/198

**Full Changelog**: https://github.com/dendrograms/astrodendro/compare/v0.2.0...v0.3.0

## 0.2.0 (2016-09-29)

- Make sure that calling structure_at with an array, list, or tuple all behave the same. [#98]
  
- Added support for linked scatter plots and multiple selections. [#104, #105, #109, #136]
  
- Added support for custom functions to define what a 'neighbor' is. [#108]
  
- Fixed a bug that caused the interactive viewer when showing a dendrogram loaded from a file. [#106, #110]
  
- Added a 'prune' method to prune dendrograms after computing them. [#111]
  
- Added support for brightness temperatures in Kelvin. [#112]
  
- Cache/memoize catalog statistics. [#115]
  
- Make sure that periodic boundaries (e.g. longitude) are properly supported. [#121]
  
- Added progress bar for catalog computation. [#127]
  
- Better support for image WCS. [#126, #137]
  
- Improve the performance of dendrogram loading. [#131]
  
- Include dendrogram parameters in HDF5 files. [#142, #145]
  
- Give HDUs names in FITS output. [#144]
  

## 0.1.0 (2013-11-09)

Initial release
