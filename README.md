# paims_codes

### Overview
This repository contains an extension to the Asymptotically Independent Markov Sampling (AIMS) algorithm for stochastic optimisation purposes and has been currently developed for the estimation of the hyper-parameters of a Gaussian process emulator. This work is part of my PhD research project at the [University of Liverpool](https://www.liv.ac.uk/risk-and-uncertainty/). Currently the only version available is on `Matlab` but it is expected to be written in `R` and `Python` in the near future. 

### Contents
* A *startup* file to load all directories needed for the sampler.  
* The `matlab-files` directory contains all files needed for the examples to run. Being  
```Matlab
    [ ... ] = parallel_aims_opt( ... ) 
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the main code for the sampler.  
* The `examples` directory contains several examples used in the paper submitted with the results obtained from such extension.  

### Citing
We would be grateful if any results based on this parallel adaptive extension are acknowledged by citing our paper. It is currently available [here](http://www.sciencedirect.com/science/article/pii/S0167947316301311)

```TeX
@Article{Garbuno2016a,
  author =   {A. Garbuno-Inigo and F.A. DiazDelaO and K.M. Zuev},
  title =    {Gaussian process hyper-parameter estimation using Parallel Asymptotically Independent Markov Sampling },
  journal =  {Computational Statistics \& Data Analysis },
  year =     {2016},
  volume =   {103},
  pages =    {367 - 383},
  doi =      {http://dx.doi.org/10.1016/j.csda.2016.05.019},
  issn =     {0167-9473},
  keywords = {Gaussian process},
  url =      {http://www.sciencedirect.com/science/article/pii/S0167947316301311}
}
```
