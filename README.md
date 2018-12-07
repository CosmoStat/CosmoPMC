# CosmoPMC
Cosmology sampling with Population Monte Carlo (PMC)

## Information

### Description

CosmoPMC is a Monte-Carlo sampling method to explore the likelihood of various
cosmological probes. The sampling engine is implemented with the package
pmclib. It is called Population MonteCarlo (PMC), which is a novel technique to
sample from the posterior (Cappé et al. 2008). PMC is an adaptive importance
sampling method which iteratively improves the proposal to approximate the
posterior. This code has been introduced, tested and applied to various
cosmology data sets in Wraith, Kilbinger, Benabed et al. (2009). Results on the
Bayesian evidence using PMC are discussed in Kilbinger, Wraith, Benabed et al.
(2010). 


### Authors

Martin Kilbinger

Karim Benabed, Olivier Cappé, Jean Coupon, Jean-François Cardoso, Gersende Fort, Henry Joy McCracken, Simon Prunet, Christian P. Robert, Darren Wraith 

### Version

1.3

### Installation

`CosmoPMC` requires the libraries `nicaea` and `pmclib`. First, download and install those packages, from their respective github pages for [nicaea](https://github.com/CosmoStat/nicaea) and [pmclib](https://github.com/martinkilbinger/pmclib).

Next, download the `CosmoPMC` package from the github repository:

```bash
git clone https://github.com/martinkilbinger/CosmoPMC
```

A new directory `CosmoPMC` will be created automatically. Change into that directory, and configure the code with the (poor-man's) python configuration script.

```bash
cd CosmoPMC
./configure.py
```

You will need to indicate paths to libraries and other flags. Type `./configure.py -h` to see all options.

After configuration, compile the code as follows:

```bash
make
```


### References

If you use CosmoPMC in a publication, please cite the last paper in the list below (Wraith, Kilbinger, Benabed et al. 2009).

[Kilbinger et al. (2011)](https://arxiv.org/abs/1101.0950): Cosmo Population Monte Carlo - User's manual. Note that earlier version of CosmoPMC <=1.2) contain `pmclib` and `nicaea` as built-in code instead of external libraries.

[Kilbinger, Benabed et al. (2012)](http://ascl.net/1212.006): ASCL link of the software package

[Kilbinger, Wraith, Benabed et al. (2010)](https://arxiv.org/abs/0912.1614): Bayesian evidence

[Wraith, Kilbinger, Benabed et al. (2009)](https://arxiv.org/abs/0903.0837): Comparison of PMC and MCMC, parameter estimation. The first paper to use CosmoPMC.

