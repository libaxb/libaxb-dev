# libaxb - Single Interface to well-tuned GPU Libraries

*libaxb* provides GPU-accelerated and vectorized linear algebra routines tuned by different vendors and developers under a single c-library interface. If you want your software to equally benefit from *cuBLAS*, *cuSparse*, *clBLAS*, *clSPARSE*, *hipSPARSE* and the many other libraries out there, *libaxb* is for you!

## Features
* C-interface for easy integration into user projects written in their favorite language
* Runtime switches for 
  -  different precision (single, double)
  -  real or complex arithmetic
  -  dense and sparse matrix types
  -  compute backends in use (e.g. switch between *cuSparse* and *clSparse* seamlessly)
* Wrappers for
  - NVIDIA's *cuBLAS* and *cuSparse*
  - AMD's *clBLAS*, *clSPARSE*, and *hipSPARSE*
  - possibly many others. Contributions welcome!
* libaxb natively provides the fast kernels from [ViennaCL](http://viennacl.sourceforge.net/)


## LICENSE
*libaxb* is available under a permissive MIT license. 
