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
