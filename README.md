# PLLPY
This is a wrapper for the [phylogenetic likelihood library](http://www.libpll.org/) (libpll) of the [Exelixis lab](http://sco.h-its.org/exelixis/index.html).

## Installation
Requirements:

- [cython](http://cython.org/)
- libpll installed from [git](https://www.assembla.com/code/phylogenetic-likelihood-library/git/nodes/master) - version 1.0.0 is too old

Installing libpll:

- on Mac PLL can be installed via [homebrew](http://brew.sh/) using the homebrew-science tap (this should also work for linuxbrew)

   - `brew tap homebrew/homebrew-science`

   - `brew install --HEAD libpll`
   
- GCC compiler versions 5 and above may need the -fgnu89-inline flag to compile libpll's inline functions, which are written in a pre-C99 style.

## Quick start:

`pllpy` has a single class, `pll`. This wraps the pllInstance type of libpll. Initialising a pll object can be done as follows:
    
    import pllpy
    instance = pllpy.pll(alignment_file, partitions, starting_tree, nthreads, random_seed)
    
- `alignment_file` is a Phylip-formatted input alignment file, passed as a string
- `partitions` is the RAxML-format partition list that associates a phylogenetic model (or models) to the alignment. This format is described in the RAxML [manual](http://sco.h-its.org/exelixis/resource/download/NewManual.pdf) - it's the format of the file passed by `-q` to RAxML. There is an example below. The value must be either a string filepath pointing to a file containing the partition description, or a string containing the partition description directly.
- `starting_tree` is one of three choices. If it is a string, then it must be a file path pointing to a newick-formatted tree file. Otherwise it should be `True`, to generate a parsimony tree, or `False`, to generate a random tree.
- `nthreads` sets the number of threads the pllInstance can use, if `pllpy` is linked to a multithreaded build of libpll.
- `random_seed` is passed to the libpll random number generator.

Maximum likelihood optimisation can be done using the `pll.optimise()` method:
    
    # Optimise tree topology and all parameters, updating the topology every 5 steps (default):
    instance.optimise(rates=True, branches=True, freqs=True, alphas=True,
                      topology=True, tree_search_interval=5, final_tree_search=True)
                      
    # Equivalently:
    instance.optimise()
                      
    # Just optimise the model parameters, not the tree topology:
    instance.optimise(topology=False)
    
    # Print the log-likelihood
    print(instance.get_likelihood())
    
## Partition format

The partition format is a list of regions of the alignment file:

    Model, label1 = Start - End
    Model, label2 = Start - End
    Model, label3 = Start - End
    
`Start - End` defines the range of alignment columns. `Model` associates a phylogenetic model to the column range. This is `GTR` for DNA data, and `DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, MTART, MTZOA, PMB, HIVB, HIVW, JTTDCMUT, FLU, STMTREV, LG4M, LG4X, GTR` for protein data. `AUTO` can be used to automatically select the best model by AIC. Adding the suffix `F` to a model uses fixed empirical base frequencies. Adding `X` lets the base frequencies be optimised.
