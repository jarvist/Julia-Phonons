# POSCAR_volumes.jl 
# Demonstration / scratchpad for using JuliaPhonons module

push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using JuliaPhonons 

for filename in ARGS
    P=read_POSCAR(open(filename))
    @printf("%s Vol: %f A^3 CubicVector: %f A\n",filename,P.volume,P.volume^(1/3))
end
