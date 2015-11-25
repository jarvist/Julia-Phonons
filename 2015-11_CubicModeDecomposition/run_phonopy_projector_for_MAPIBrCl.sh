for i in MAPbBr MAPbCl MAPbI
do
 cd "${i}"
 julia ../phonopy_projector.jl > "../${i}.dat"
 cd -
done
