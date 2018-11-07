# MAPI_acoustic_BZ.jl

push!(LOAD_PATH,"../../src") # Temporary versions of modules in PWD
using JuliaPhonons 

P=read_POSCAR(open("../POSCAR"),expansion=[4,4,4])
eigenvectors,eigenmodes=read_meshyaml(open("../mesh.yaml"),P)

output_conflated_xyz(P,0,(eigenvectors[1],eigenvectors[1]), 
                    (1,1.0), 
                    steps=128, repeats=1, q=[0.25,0.25,0.25]) # generates files anim_{count}.xyz

output_conflated_xyz(P,1,(eigenvectors[2],eigenvectors[2]), 
                    (1,1.0), 
                    steps=128, repeats=1, q=[0.25,0.25,0.25]) # generates files anim_{count}.xyz

output_conflated_xyz(P,2,(eigenvectors[3],eigenvectors[3]), 
                    (1,1.0), 
                    steps=128, repeats=1, q=[0.25,0.25,0.25]) # generates files anim_{count}.xyz


end
output_conflated_xyz(P,0,(eigenvectors[1], eigenvectors[2],eigenvectors[3]), 
                    (1,1,1), 
                    steps=64, repeats=1, q=[0.25,0.25,0.25]) # generates files anim_{count}.xyz


