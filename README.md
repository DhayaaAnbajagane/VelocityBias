# Velocity Dispersion/Bias Scaling Relations from Cosmological Simulations

This repository hosts both the data and corresponding interpolation scripts for the galaxy/DM velocity dispersion and galaxy velocity bias scaling relation parameters of halos in five simulations IllustrisTNG, Magneticum Pathfinder, Bahamas + Macsis, The Three Hundred, and MultiDark Planck-2 from redshifts 0 < z < 1. We also provide a simulation ensemble-estimated theoretical prior (mean and 1sigma uncertainty) on the velocity bias. The scaling parameters are computed as described in [Anbajagane et. al 2021](https://arxiv.org/abs/2110.01683). **Citation to this publication is requested if these scaling parameters are used.**

The scaling relations with halo mass are available for the following properties:

1. sigma_DM: The DM velocity dispersion
2. sigma_Sat: The satellite galaxy velocity dispersion
3. b_v: The velocity bias


## QUICKSTART

```

import sys
sys.path.append("<path to 'VelocityBias'>/VelocityBias")

import Interpolator
import numpy as np

print("List of available sims is " + str(Interpolator.avail_sims))
print("List of available scaling parameters is " + str(Interpolator.avail_params))

M200c = 10**np.linspace(13.5, 14.5, 100)

sigmaDM  = Interpolator.sigmaDM(M200c, z = 0.8, parameter = 'mean', sim = 'TNG')
                                     
sigmaSat = Interpolator.sigmaSat(M200c, Mstarsat_Th = 1e10, z = 0.3, parameter = 'mean', sim = 'MGTM')
                                
bv = Interpolator.velocity_bias(M200c, Mstarsat_Th=1e10, z=0.2, sims = ['BM', 'TNG', 'MGTM', 'The300'])

```

### Caveats

1. MultiDark Planck-2 estimates are only available for sigmaSat and not for sigmaDM or b_v.
2. The slope and scatter output from the sigmaSat function are point estimates rather than mass-dependent ones. So the output at a given stellar mass threshold and redshift will be the same regardless of the input mass.
3. For the velocity bias, only the mean relation is available. While the slope can be estimated from the mean relation, the code does not compute it. The scatter about this relation was not estimated.

If you find any errors/bugs in the code, please reach out to dhayaa@uchicago.edu
