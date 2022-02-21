# FailureOfInhibition2022

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> FailureOfInhibition2022

To (locally) reproduce this project, do the following:

1. Download this code base (preferably via `git clone`)
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> pkg"registry add https://github.com/grahamas/SmithThesisRegistry.git"
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box. 

*NOTE*: this adds a new registry containing packages written expressly for this paper and related thesis.
