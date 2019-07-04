# MAMBO_Functions

Codes for the functions that were used in the MAMBO algorithm described by Garza et al. (2018). These functions were written on Cython and are intended to be used integrated with COBRApy (https://github.com/opencobra/cobrapy) on a python2.7 platform, using a collection of genome-scale metabolic models with model SEED identifiers for metabolites and reactions.  

The script "bottom_up_ecology.pyx" contains the functions for: 1) identifying exchange reactions, 2) setting-up an initial environment 3) Performing a Markov-Chain step to determine the accept/reject probabilities for random changes to the environment. In order to use these functions, one must install the prerequisites listed below and compile the script using the "setup_bte_functions.py" script.

The functions are written on cython, thus they require compilation. To do this, one can use the terminal command to compile to their working directory:

$ setup.py build_ext --inplace

Other ways of compiling Cython code should also be available. [Check here](http://cython.readthedocs.io/en/latest/src/reference/compilation.html)

After compilation, the script may be imported into a python script by "import bottom_up_ecology". An example of usage for maximizing the correlation between steady-state solutions for biomass functions and the relative abundance of species determined by metagenomic experiments is given below. For more details, see Garza et al.(2016).

###Example of MCMC optimization
for this example, a python environment is assumed. It is also assumed that one has a group of GSMMs on a python list named "models" and, in the same order, a numpy array of numbers with the relative abundance of bacteria named "relative_abundance"  


import bottom_up_ecology_functions as bte


evolution_of_exchange_reactions = {} # A python dictionary to store the results for the MCMC

initial_environment = bte.starting_environment(models) # defines an initial random environment

metabs = initial_environment.keys() # defines the metabolites that will conform the search-space.

for i in xrange(1000): # 1000 Markov chain optimization steps

    bte.single_MCMC_run(numpy.random.random_integers(0, len(metabs)), models, initial_environment, metabs, relative_abundance, evolution_of_exchange_reactions , i)


### Dependencies

List of dependencies with the versions that were used for the data described on the manuscript:

- COBRApy (v. 0.4.1)
- Numpy (v. 1.11.0)
- Gurobi (v. 6.5)
- Cython (v. 0.24)
- Scipy (v. 0.17.1)


### Models

1562 genome-scale metabolic models were used in the study. These models can be obtained from this repository, by downloding the files:

modelsaa
md5: 8f0d73a450048b2a8aab6e7b9f985b4d

modelsab
md5: e386510345c950ace6ec66a85c488ebb  

modelsac
md5: c4d2ed40d4a77de19719c9b19fa9e1e4

modelsad
md5: 7f204c1c41b6fa4013e9aab63a2756fd

modelsaf
md5: bd7c0c3ebc8383a1d3653d5d2bd468cf

once these files are obtained, the linux command:

$ cat models*>models.tar.gz

will join them into a single file.

models.tar.gz
md5: b8694cc4da915317a5549ab46f0a3032

