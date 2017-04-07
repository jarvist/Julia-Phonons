# set grid_mode, 1

alter elem I, vdw=1.0
alter elem Pb, vdw=1.0
show spheres, elem I
show spheres, elem Pb

rebuild spheres

# redefine 'hydrogen bond' chunky dashed yellow lines into something nicer
set dash_gap, 0.0
set dash_radius,0.02
# draw distance-labelled bonds between Pb and Is
dist name Pb*, name I*,mode=3,cutoff=4
dist name I*,name I*, mode=3, cutoff=5

# Ray trace defaults
bg_color white
show spheres
set ray_trace_fog,0
set ray_shadows,0
set antialias,1
set spec_reflect, 0.0
set reflect, 0.5


