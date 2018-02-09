# defmod-swpc -- A finite element and finite difference mixed code for multiphysics and multiscale crustal motions. 

Corresponding Author: Chunfang Meng ( cmeng (at) mit.edu )   
Contributors: Tabrez Ali ( tabrez.ali (at) gmail.com ) and Takuto Maeda ( maeda (at) eri.u-tokyo.ac.jp )

* * *

## DESCRIPTION
This code combines two open source codes Defmod(-dev) and OpenSWPC, forming an FE-FD mixed code. 

To install and use the "pure" FE or FD code, visit   
https://bitbucket.org/stali/defmod    
https://github.com/takuto-maeda/OpenSWPC

To install and run this "mixed" code, see doc/INSTALL

## FEATURES
1. Forward incremental Lagrange Multiplier method for modeling fault slip honoring slip weakening and rate-and-state frictional laws, verified with SCEC benchmark problems.
2. Adaptive (Quasi-)static-dynamic "hybrid" (FE) solver that covers multi-temporal scales, e.g. fault slip (fast) triggered by pore pressure changes (slow).  
3. Finite element and finite difference coupling that covers multi-spatial scales, e.g. localized fault slip and resulting ground motions. 
4. Generalized-alpha (implicit dynamic) solver that is unconditionally stable, i.e. does not suffer CFL restriction, while allowing fault rupture.  

* * *

## LICENSE
MIT License, see LICENSE for details.

The authors appreciate that the users cite the following papers in any publications employing this code. For feedback or contribution, please contact the corresponding author. 


* * *

## REFERENCES

Chunfang Meng, Benchmarking Defmod, an open source FEM code for modeling episodic fault rupture, Computers & Geosciences, Volume 100, March 2017, Pages 10-26, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2016.11.014.

Takuto Maeda, Shunsuke Takemura and Takashi Furumura (2017), OpenSWPC: An open-source integrated parallel simulation code for modeling seismic wave propagation in 3D heterogeneous viscoelastic media, Earth, Planets and Space, 69:102, https://doi.org/10.1186/s40623-017-0687-2. 

C. Meng, H. Wang, A finite element and finite difference mixed approach for modeling fault rupture and ground motion, Computers & Geosciences 113 (2018) 54 – 69, https://doi.org/10.1016/j.cageo.2018.01.015.
