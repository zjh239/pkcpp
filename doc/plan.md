
### RDF analysis

1. From binning the atoms, to dealing with the periodic bins, the code from neighbor finder can be reused.

2. The different part is that we do not construct the neighbor list, but put the pair distance into a histogram according to pair type.

First, the pair distance data is calculated using the bin mesh function. The data is put into a list according to the atom type of the pair.

So that, the counting list need to be initialized first based on the number of atom types.
