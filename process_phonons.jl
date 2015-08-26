import YAML

mesh = YAML.load(open("mesh.yaml"))
POSCAR = YAML.load(open("POSCAR.yaml"))

NATOMS=12 # Should _not_ be hardcoded!

println(POSCAR["positions"])
positions=[POSCAR["positions"][n][d] for n=1:NATOMS,d=1:3 ]

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
# OK; we've built atomnames[] to reproduce POSCAR

mode=1
# Data structure looks like: mesh["phonon"][1]["band"][2]["eigenvector"][1][2][1]
for (eigenvector,freq) in mesh["phonon"][1]["band"]
    filename= @sprintf("anim_%02d.xyz",mode)
    anim=open(filename,"w")
    println("freq ==> ",freq[2])
#    println("phonon[\"eigenvector\"] ==>",phonon["eigenvector"])
#    println("eigenvector ==> ",eigenvector[2])
#    for atom in eigenvector[2]
#        #println ("atom ==> ",atom)
#        @printf("atom %e %e %e \n",atom[1][1],atom[2][1],atom[3][1])
#        disp=[atom[1][1],atom[2][1],atom[3][1]] #pull data out of YAML into sensible Julia Vector
#        println(disp)
#    end

    realeigenvector=[ eigenvector[2][n][d][1] for n=1:NATOMS, d=1:3 ]
    # Array comprehension to reform mesh.yaml format into [n][d] shap
    realeigenvector=reshape(realeigenvector,NATOMS,3) # doesn't do anything?

    for phi=0:pi/8:2*pi-1e-6 #slightly offset from 2pi so we don't repeat 0=2pi frame
        projection=6.3*positions+10*realeigenvector*sin(phi) # this does all the maths

        # output routines to .xyz format
        @printf(anim,"%d\n\n",NATOMS) # header for .xyz multi part files
        for i=1:size(projection,1) 
#            println("u ==> ",projection[i,:])
            @printf(anim,"%s %f %f %f\n",atomnames[i],projection[i,1],projection[i,2],projection[i,3])
        end
    end
    close(anim)
    mode=mode+1
end

    #POSCAR["positions"][1]
