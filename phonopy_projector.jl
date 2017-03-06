# phonopy_project.jl
# Demonstration / scratchpad for using JuliaPhonons module

push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using JuliaPhonons 

P=read_POSCAR(open("POSCAR"),expansion=[2,2,2])
eigenvectors,eigenmodes=read_meshyaml(open("mesh.yaml"),P)

myf=open("mode_decomposition.dat","w")
gnuplot_header(P,f=myf)

for (count,(eigenvector,eigenmode)) in enumerate(zip(eigenvectors,eigenmodes))
    output_animated_xyz(P,count,eigenvector,eigenmode) # generates files anim_{count}.xyz
    
    decompose_eigenmode_atomtype(P,count,eigenvector,eigenmode)
    # Output again, to a file for later printing
    decompose_eigenmode_atomtype(P,count,eigenvector,eigenmode,f=myf)

    decompose_eigenmode_atom_contributions(P,count,eigenvector)
end

close(myf)
