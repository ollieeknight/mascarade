# mascarade (development version)
* Modified `generateMask()` to ensure only one label appears per cluster, even when a cluster has multiple disconnected regions. The largest region of each cluster receives the label.
* Added `label` column to the output of `generateMask()` to control which parts get labeled.

# mascarade 0.3.1
* fancyMask support setting colors manually or to inherit from the main layer.

# mascarade 0.3.0

* Significant performance improvements
* Introduced `generateMaskSeurat()` and `fancyMask()` helper functions.

# mascarade 0.2.0

* Major rewrite of generateMask function, with a change of interface
