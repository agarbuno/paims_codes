# paims_codes

### Overview
This repository contains an extension to the Asymptotically Independent Markov Sampling (AIMS) algorithm for stochastic optimisation purposes and has been currently developed for the estimation of the hyper-parameters of a Gaussian process emulator. This work is part of my PhD research project at the [University of Liverpool](https://www.liv.ac.uk/risk-and-uncertainty/). Currently the only version available is on `Matlab` but it is expected to be written in `R` and `Python` as well in the near future. 

### Contents
* A *startup* file to load all directories needed for the sampler.  
* The `matlab-files` directory contains all files needed for the examples to run. Being  
```Matlab
    [ ... ] = parallel_aims_opt( ... ) 
```
     the main code for the sampler. 
* The `examples` directory contains several examples used in the paper submitted with the results obtained from such extension.  

### Citing
We would be grateful if any results based on the this parallel adaptive extension are acknowledge by citing our paper. It is currently available on [arxiv](http://arxiv.org/abs/1506.08010)

```TeX
Article{
  Title                    = {{Gaussian process hyper-parameter estimation using parallel asymptotically
independent Markov sampling}},
  Author                   = {A. Garbuno-Inigo, F.A. DiazDelaO, K.M. Zuev},
  Journal                  = {Arxiv},
  Year                     = {2015},
  Archiveprefix            = {arXiv},
  Arxivid                  = {1506.08010},
  Url                      = {http://arxiv.org/abs/1506.08010}
}
```
