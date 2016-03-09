# phonopy_project.jl
# Demonstration / scratchpad for using JuliaPhonons module

push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using JuliaPhonons 

P=read_POSCAR(open("POSCAR"))
eigenvectors,eigenmodes=read_meshyaml(open("mesh.yaml"),P)

for (count,(eigenvector,eigenmode)) in enumerate(zip(eigenvectors,eigenmodes))
    output_animated_xyz(P,count,eigenvector,eigenmode) # generates files anim_{count}.xyz
    decompose_eigenmode_atomtype(P,count,eigenvector,eigenmode)
    decompose_eigenmode_atom_contributions(P,count,eigenvector)
end
