## Circuitscape processing..
library(JuliaConnectoR)

## read in standard Circuitscape.ini file
# use grep and gsub to replace specific text strings
# need to also use taxonpath
# line 9: habitat_file = taxonpath/resistance.tif
## or substitute resistance file derived from HDM (or HIM)
# line 18 polygon_file = orphans.tif
# line 30 mask_file = ??
# line 39 point_file = taxonpath/preclusters.tif
# line 47 output_file = taxonpath

