# phonopy_project.jl
# Demonstration / scratchpad for using JuliaPhonons module

push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using JuliaPhonons 

P=read_POSCAR(open("POSCAR"),expansion=[2,2,2])
eigenvectors,eigenmodes=read_meshyaml(open("mesh.yaml"),P)

myf=open("mode_decomposition.dat","w")
gnuplot_header(P,f=myf)

output_conflated_xyz(P,1,(eigenvectors[4],eigenvectors[12],eigenvectors[18]),(eigenmodes[4],eigenmodes[12],eigenmodes[18]) ) # generates files anim_{count}.xyz

output_conflated_xyz(P,2,(eigenvectors[4],eigenvectors[7],eigenvectors[10],eigenvectors[13],eigenvectors[16],eigenvectors[19]),(eigenmodes[4],eigenmodes[7],eigenmodes[10],eigenmodes[13],eigenmodes[16],eigenmodes[19]) ) # generates files anim_{count}.xyz
 
output_conflated_xyz(P,3,(eigenvectors[4],eigenvectors[7],eigenvectors[10],eigenvectors[13],eigenvectors[16],eigenvectors[19],eigenvectors[22],eigenvectors[25],eigenvectors[28],eigenvectors[31],eigenvectors[34]),(eigenmodes[4],eigenmodes[7],eigenmodes[10],eigenmodes[13],eigenmodes[16],eigenmodes[19],eigenmodes[22],eigenmodes[25],eigenmodes[28],eigenmodes[31],eigenmodes[34]) ) # generates files anim_{count}.xyz
 
for (count,(eigenvector,eigenmode)) in enumerate(zip(eigenvectors,eigenmodes))
    output_animated_xyz(P,count,eigenvector,eigenmode) # generates files anim_{count}.xyz
 
    decompose_eigenmode_atomtype(P,count,eigenvector,eigenmode)
    # Output again, to a file for later printing
    decompose_eigenmode_atomtype(P,count,eigenvector,eigenmode,f=myf)

    decompose_eigenmode_atom_contributions(P,count,eigenvector)
end

close(myf)
