#!/bin/julia
# What Would Hamilton Do? (WWHD?)
using Quaternions

# Origin of MA in Fract Coords
origin=[0.5 0.5 0.5]' # origin of MA in unit cell, fractional coords

POSCARlines = readlines(open("POSCAR"))

shift=[0.05*randn() 0.05*randn() 0.05*randn()]'

# Generate a normalized Quaternion with normally distributed random 4 component
qr=normalize(Quaternion(randn(),randn(),randn(),randn()))
# This should be an evenly distributed rotation matrix
rotate=rotationmatrix(qr)

#rotate=eye(3) #3-dimensional identity matrix; null rotation operation

#println("Determinant (should be 1): ",det(rotate))
#println("Does Q^TQ = I? \n",rotate*rotate')

for POSCARline in enumerate(POSCARlines)
    # Start line 10 if Selective dynamics; 9 otherwise
    if (POSCARline[1]>=9 && POSCARline[1]<=9+8) #if MA; for Cubic Perov structs.
        r=float(split(POSCARline[2])) # "1.2 3.4 5.6" --> [1.2; 3.4; 5.6]
        r=shift+origin+(rotate*(r-origin)) # apply shift + rotate, shift back
        r=mod.(r,1) # Modulo arithmatic so that all {x,y,z} values are on [0,1]
        @printf "  %.16f  %.16f  %.16f \n" r[1] r[2] r[3]
        # 16 digits of precision seems to match VASP output
    else # otherwise echo line unchanged to output
        println(POSCARline[2])
    end
end

