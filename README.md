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
5. (Partially) synchrize FE domain pressure/permeability with a finite volume code [pflotran](https://bitbucket.org/pflotran/pflotran/wiki/Home).

## Example descriptions
* F3D, slip on curved fault induced by increasing differential loading, 3D domain.
* F2Dp, production fluctuation and rate-and-state friction scenarios, 2D domain.
* F3Dp, production induced rupture on curved fault, 3D domain.
* F3Db (F3Db\_usg for unstructured and shared FV/FE mesh), minimalistic example of injection induced fault slip to demonstrate FV(pflotran)->FE(defmod)->FD(swpc) binding, 3D domain, transient and tensor valued permeability. To activate FV-FE binding, launch defmod with -fv [argv]:
    * -fv 1: structured FV single phase;
    * -fv 2: structured FV multiphase;
    * -fv 3: unstructured FV single phase;
    * -fv 4: unstructured FV multiphase.
* HF3D, injection induced rupture on splay fault, 3D domain (modified from SCEC14/15).
* FE, legacy examples including Mandel (2D/3D) benchmarks.
* SCEC, SCEC problems number 10, 102 and 205, waveforms produced by both the FE and FD modules.
* F3DX, Modified SCEC problems number 14 and 15 to have intersecting faults.
* SCEC2D, SCEC 2D problems number 10, 11 and 102 to compare hybrid (implicit-explicit) models against generalized-alpha (implicit dynamic) models. 
* CR3D, a 3D equivalent to [Cappa and Rutqvist 2011](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2011GL048487) model, demonstrating the FV->FE->FD workflow, follow [README](https://github.com/Chunfang/defmod-swpc/blob/master/example/CR3D) to execute.\
Examples F3D, F3Db, F2Dp, SCEC2D and FE can run on a desktop computer, while others shall run on a cluster. The execution commands are equivalent for both the desktop and cluster, see <model>.sh files for details.

* * *

## LICENSE
MIT License, see LICENSE for details.

The authors appreciate that the users cite the following papers in any publications employing this code. For feedback or contribution, please contact the corresponding author. 


* * *

## REFERENCES
Chunfang Meng, Benchmarking Defmod, an open source FEM code for modeling episodic fault rupture, Computers & Geosciences, Volume 100, March 2017, Pages 10-26, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2016.11.014.

Takuto Maeda, Shunsuke Takemura and Takashi Furumura (2017), OpenSWPC: An open-source integrated parallel simulation code for modeling seismic wave propagation in 3D heterogeneous viscoelastic media, Earth, Planets and Space, 69:102, https://doi.org/10.1186/s40623-017-0687-2. 

C. Meng, H. Wang, A finite element and finite difference mixed approach for modeling fault rupture and ground motion, Computers & Geosciences 113 (2018) 54 – 69, https://doi.org/10.1016/j.cageo.2018.01.015.
