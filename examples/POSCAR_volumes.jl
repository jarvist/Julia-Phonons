# POSCAR_volumes.jl 
# Demonstration / scratchpad for using JuliaPhonons module

push!(LOAD_PATH,"../src") # Temporary versions of modules in PWD
using JuliaPhonons 

# OK, a pretty weird thing to do!
# But the POSCAR reader, part of JuliaPhonons, can also do basic stats on the unit cell.

for filename in ARGS
    P=read_POSCAR(open(filename))
    @printf("%s Vol: %f A^3 CubicVector: %f A\n",filename,P.volume,P.volume^(1/3))
end
