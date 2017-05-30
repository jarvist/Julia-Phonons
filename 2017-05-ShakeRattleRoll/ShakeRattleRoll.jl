# phonopy_project.jl
# Demonstration / scratchpad for using JuliaPhonons module

push!(LOAD_PATH,"../src") # Temporary versions of modules in PWD
using JuliaPhonons 

# The REPO comes with a basic MAPI cubic unit cell Phonopy calculation
P=read_POSCAR(open("../2015-11_CubicModeDecomposition/MAPbI/POSCAR"),expansion=[2,2,2])
eigenvectors,eigenmodes=read_meshyaml(open("../2015-11_CubicModeDecomposition/MAPbI/mesh.yaml"),P)

myf=open("mode_decomposition.dat","w")
gnuplot_header(P,f=myf)

mode=16

println("My eigenvector $mode, ")
show(eigenvectors[mode])
println("My eigenvalue for $mode, ")
show(eigenmodes[mode])
println()

# See: https://stackoverflow.com/a/11319919
# Essentially gen original + perturbed structure; best Quaternions minimises error on forwards / backwards transition
#

# Eeek! Hard coded.
println("C: ",eigenvectors[mode][1,:])
println("N: ",eigenvectors[mode][2,:])

end

for (count,(eigenvector,eigenmode)) in enumerate(zip(eigenvectors,eigenmodes))
    output_animated_xyz(P,count,eigenvector,eigenmode) # generates files anim_{count}.xyz
 
    decompose_eigenmode_atomtype(P,count,eigenvector,eigenmode)
    # Output again, to a file for later printing
    decompose_eigenmode_atomtype(P,count,eigenvector,eigenmode,f=myf)

    decompose_eigenmode_atom_contributions(P,count,eigenvector)
end

close(myf)
