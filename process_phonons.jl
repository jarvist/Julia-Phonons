import YAML

# This would be better in a general library!
atomicmass = Dict{AbstractString,Float64}("H"=>1, "C" => 12.01, "N" => 14.01, "S" => 32.07, "Zn" => 65.38, "I" => 126.9, "Pb" => 207.2)
#Zn is actually Sn; stupid work-around for Pymol seeing 'Sn' as 'S'
#The initial test cases which I am covering, are phonons of CH3NH3.PbI3 and SnS

mesh = YAML.load(open("mesh.yaml"))     #Phonopy mesh.yaml file; with phonons
POSCAR = YAML.load(open("POSCAR.yaml")) #A bit of a twisted POSCAR->yaml format

NATOMS=sum(POSCAR["speciescount"]) # Reads atoms from sum of POSCAR species line 
println("NATOMS ==> ",NATOMS)

positions=[POSCAR["positions"][n][d]::Float64 for n=1:NATOMS,d=1:3 ] # Array comphrehension to eval as floats
println("positions ==> ",positions)

lattice=[ POSCAR["lattice"][d][a]::Float64 for d=1:3,a=1:3] # Nb: Array comp for sensible lattice[a][b] julia type 
println("lattice ==> ",lattice)

# SUPERCELL definition
supercellexpansions=[ a*lattice[1,:] + b*lattice[2,:] + c*lattice[3,:] for a=0:1,b=0:1,c=0:1 ] #generates set of lattice vectors to apply for supercell expansions
println("supercellexpansions ==>",supercellexpansions)

# The following is probably overkill, but reads the VASP atom formats + expands
# into a one dimensional string vector 
#     species:    C    N    H    Pb   I
#     speciescount:   1     1     6     1     3
atomnames=String[]
println(POSCAR["speciescount"])
println(POSCAR["species"])

for (speciescount,species) in zip(POSCAR["speciescount"],POSCAR["species"])
    println("speciescount => ",speciescount," species => ",species)
    for i=1:speciescount push!(atomnames,species) end
    println(atomnames)
end
# OK; we've built atomnames[], mainly for use with outputs... 

function output_animated_xyz(eigenmode,eigenvector,freq,steps=32)
    filename= @sprintf("anim_%02d.xyz",eigenmode)
    anim=open(filename,"w")

    for phi=0:2*pi/steps:2*pi-1e-6 #slightly offset from 2pi so we don't repeat 0=2pi frame
#        projection= lattice[1][1]*positions + 2*realeigenvector*sin(phi) # this does all the maths
#        println("Lattice: ",lattice,"\n Eigenvec: ",realeigenvector)
        #projection=positions*lattice + 2*sin(phi)*realeigenvector
#        println("Projection: ",projection)

        # output routines to .xyz format
        @printf(anim,"%d\n\n",NATOMS*length(supercellexpansions)) # header for .xyz multi part files
        for i=1:NATOMS
#            println("u ==> ",projection[i,:])
#            println("Positions[",i,"]: ",positions[i,:])
#            println("Realeigenvector[",i,"]: ",realeigenvector[i,:])

            projection=lattice * (positions[i,:]' + 0.02 * sqrt(atomicmass[atomnames[i]]) * eigenvector[i,:]'*sin(phi))
            for supercellexpansion in supercellexpansions
                supercellprojection=projection+supercellexpansion'
                @printf(anim,"%s %f %f %f\n",atomnames[i],supercellprojection[1],supercellprojection[2],supercellprojection[3])
            end
        end
    end
    close(anim)


end

# Data structure looks like: mesh["phonon"][1]["band"][2]["eigenvector"][1][2][1]
for (eigenmode,(eigenvector,freq)) in enumerate(mesh["phonon"][1]["band"])
    println("freq (THz) ==> ",freq[2], "Wavenumbers (cm-1, 3sf) ==> ",freq[2]*33.36)
#    println("phonon[\"eigenvector\"] ==>",phonon["eigenvector"])
#    println("eigenvector ==> ",eigenvector[2])
#    for atom in eigenvector[2]
#        #println ("atom ==> ",atom)
#        @printf("atom %e %e %e \n",atom[1][1],atom[2][1],atom[3][1])
#        disp=[atom[1][1],atom[2][1],atom[3][1]] #pull data out of YAML into sensible Julia Vector
#        println(disp)
#    end

    realeigenvector=[ eigenvector[2][n][d][1]::Float64 for n=1:NATOMS, d=1:3 ]
    # Array comprehension to reform mesh.yaml format into [n][d] shap
    realeigenvector=reshape(realeigenvector,NATOMS,3) # doesn't do anything?

    output_animated_xyz(eigenmode,realeigenvector,freq)

    normsum=0.0
    normsummassweighted=0.0
    for i=1:NATOMS
        normsum+=norm(realeigenvector[i,:])
        normsummassweighted+=norm(realeigenvector[i,:])/sqrt(atomicmass[atomnames[i]])
    end
    println("Norm sum: ",normsum, " Norm sum(mass weighted): ",normsummassweighted)
    for i=1:NATOMS    
        println("Mode: ",eigenmode," Atom: ",i," ",atomnames[i],
        "\n",
        " Norm: ", norm(realeigenvector[i,:]),
        " Norm(weighted): ",norm(realeigenvector[i,:])/normsum,
        " Norm(mass weighted): ",(norm(realeigenvector[i,:])/sqrt(atomicmass[atomnames[i]]))/normsummassweighted)
    end
    
#=    
    for I=10:12 # iodine indexes, hard coded
        println(show(positions))
        println(show(positions[I,:]))
        PbI=positions[I,:]-positions[9,:] # Vector from Pb to I
        show(PbI)
        @printf("I: %d PbI: %f %f %f Phonon-u: %f %f %f PbI.Phonon-u: %f\n",I,
        PbI[1],PbI[2],PbI[3],
        realeigenvector[I,1],realeigenvector[I,2],realeigenvector[I,3],
        dot(PbI::Float64,realeigenvector[I]::Float64) )
    end
=#
end

