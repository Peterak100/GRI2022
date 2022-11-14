## Circuitscape processing..
# library(JuliaConnectoR) ? use this package to run Julia commands?

## read in standard Circuitscape.ini file
# use grep and gsub to replace specific text strings
# need to also use taxonpath
# line 9: habitat_file = taxonpath/resistance.tif
## or a substitute resistance file derived from HDM or HIM
# line 18 polygon_file = orphans.tif?
## TRIED THIS - IT WORKS BUT MAKES NEGLIGIBLE DIFFERENCE SO NOT WORTH IT
# line 30 mask_file = ??
## MASK FILE NEEDS TO be in tab-delimited text with a .txt extension
##  CHECK THAT IT IS PLAIN TEXT FORMAT WITH EXTRANEOUS "" MARKS
##  TRIED THIS - CIRCUITSCAPE CRASHED
# line 39 point_file = taxonpath/preclusters.tif
# line 47 output_file = taxonpath

## HOW TO USE CHOLMOD?? Only works with double precision??


## 1st Circuitscape run in Julia --------
## cd /home/peter/data
# julia
# using Circuitscape
## either:
# compute("/home/peter/data/taxa/taxon/Circuitscape_custom1.ini")
## or:
# for i in eachline("batch_jobs.txt")
# compute(replace("/home/peter/data/taxa/taxon987/Circuitscape_custom1.ini",
#                "taxon987" => i))
# end

## or for testing:
# for i in eachline("batch_jobs.txt")
# println(replace("/home/peter/data/taxa/taxon987/Circuitscape_custom1.ini",
#                "taxon987" => i))
# end

## 2nd Circuitscape run in Julia --------
## cd /home/peter/data
# julia
# using Circuitscape
## either:
# compute("/home/peter/data/taxa/taxon/Circuitscape_custom2.ini")
## or:
# for i in eachline("batch_jobs.txt")
# compute(replace("/home/peter/data/taxa/taxon987/Circuitscape_custom2.ini",
#                "taxon987" => i))
# end

### both all-to-one and one-to-all methods failed for the above!!

## Run from Julia singly:
# using Circuitscape
# compute("/home/peter/data/taxa/
#         Pseudemoia_cryodroma/Circuitscape_custom2.ini")
# 18 patches (153 pairwise comparisons) took just over 10 mins in Windows
#  and ~11.5 mins in Ubuntu

