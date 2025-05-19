# Process parameters
scanning_speed=0.5 # m/s
power=24 # W (this is the effective power so multiplied by eta)
R=70e-6 # m (this is the effective radius)

# Geometric parameters
thickness=100e-6 # m
xmin=-80e-6 # m
xmax=240e-6 # m
ymin=${fparse -thickness}
surfacetemp=300 # K (temperature at the other side of the plate)

# Time stepping parameters
endtime=240e-6 # s
timestep=${fparse endtime/1000} # s
