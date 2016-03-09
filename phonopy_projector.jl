# phonopy_project.jl
# Demonstration / scratchpad for using JuliaPhonons module

push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using JuliaPhonons 

P=read_POSCAR(open("POSCAR"))
read_meshyaml(open("mesh.yaml"),P)
