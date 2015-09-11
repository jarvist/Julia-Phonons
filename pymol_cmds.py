set grid_mode, 1

bond elem Pb, elem I
bond elem I, elem I

alter elem I, vdw=1.0
alter elem Pb, vdw=1.0
show spheres, elem I
show spheres, elem Pb

rebuild spheres

