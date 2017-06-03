# conflated_MAPI_example.jl
# Demonstration / scratchpad for Conflated Phonons in MAPI
#
# Used to generate this YouTube video,
# https://www.youtube.com/watch?v=7UoPgFJcRoI
#
# I've been playing with a new method to show multiple phonon modes in
# a material simultaneously. I've called it 'conflation'; essentially it tries
# to show multiple phonon modes simultaneously, with the correct relative
# frequency of motion. The result looks quite similar to molecular dynamics.
#
# The material is Methyl-ammonium lead-halide perovskite, a material of
# interest for solution processed solar cells. The underlying data comes from
# a Phonopy calculation with force-constants from a VASP plane-wave
# density-functional-theory calculations. For the 1x1x1 pseudo-cubic unit cell
# here, this generates a set of 36 eigenvectors and eigenmodes, but due to the
# (slightly broken) cubic symmetry, sets of 3 within those modes are
# approximately equivalent.
#
# This is a representation of all inorganic-cage modes, and the first
# organic-cation mode. [Modes 4,7,10,13,16,19 .]
#
# The render is done with the marvellous VMD, and the Tachyon ray tracer.
#
# The modes are shown simultaneously, with the total motion scaled relative to
# the energy (frequency) of the phonon mode. Within the background of
# 8 oscillations of the lowest frequency mode (~1 Thz), the other modes are
# sped-up relative to their actually frequency. In the first frame of the
# movie, the structure is in it's ground state, but as the eigenmodes are not
# harmonic within each other, the material will never return to this structure.
#
# In total, the movie shows about 8 picoseconds of real time.

push!(LOAD_PATH,"../src") # Temporary versions of modules in PWD
using JuliaPhonons 

P=read_POSCAR(open("POSCAR"),expansion=[4,4,2])
eigenvectors,eigenmodes=read_meshyaml(open("mesh.yaml"),P)

#output_conflated_xyz(P,1,(eigenvectors[4],eigenvectors[12],eigenvectors[18]),(eigenmodes[4],eigenmodes[12],eigenmodes[18]) ) # generates files anim_{count}.xyz

output_conflated_xyz(P,1,(eigenvectors[1],eigenvectors[2]), (eigenmodes[4],eigenmodes[4]) , bz=[0.1,0,0]) # generates files anim_{count}.xyz
 
output_conflated_xyz(P,2,(eigenvectors[4],eigenvectors[7],eigenvectors[10],eigenvectors[13],eigenvectors[16],eigenvectors[19]),(eigenmodes[4],eigenmodes[7],eigenmodes[10],eigenmodes[13],eigenmodes[16],eigenmodes[19]) , sound=true,bz=[0,0,1]) # generates files anim_{count}.xyz
 
output_conflated_xyz(P,3,
(eigenvectors[1],eigenvectors[10],eigenvectors[13],eigenvectors[16],eigenvectors[19]), 
(eigenmodes[4],eigenmodes[10],eigenmodes[13],eigenmodes[16],eigenmodes[19]) , 
bz=[1,0,0]) # generates files anim_{count}.xyz
 
#output_conflated_xyz(P,3,(eigenvectors[4],eigenvectors[7],eigenvectors[10],eigenvectors[13],eigenvectors[16],eigenvectors[19],eigenvectors[22],eigenvectors[25],eigenvectors[28],eigenvectors[31],eigenvectors[34]),(eigenmodes[4],eigenmodes[7],eigenmodes[10],eigenmodes[13],eigenmodes[16],eigenmodes[19],eigenmodes[22],eigenmodes[25],eigenmodes[28],eigenmodes[31],eigenmodes[34]) ) # generates files anim_{count}.xyz
