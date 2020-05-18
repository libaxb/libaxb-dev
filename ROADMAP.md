
# Roadmap for libaxb

## Provide operations from CUDA libraries

Currently only the fallback kernels are available for matrix and vector operations.
libaxb becomes much more attractive once the operations from the CUDA libraries (CUBLAS, CUSPARSE) are hooked in.

## Provide operations from AMD libraries

Same rationale as before.

## Provide operations from Intel GPU libraries

Same rationale as before.

## Integration into PETSc

In order to quickly get more users of libaxb, there is a [libaxb pull request for PETSc](https://gitlab.com/petsc/petsc/-/merge_requests/2073) in progress.
Finish the pull request.

## Provide a file-based kernel selector

Users should be able to specify which kernels to use from text files (or options).
This could be a JSON format.

## Benchmarks

Demonstrate the usefulness of the mix-and-match approach from PETSc.
In particular, show how easy it is to mix different CUDA kernels from different libraries.

