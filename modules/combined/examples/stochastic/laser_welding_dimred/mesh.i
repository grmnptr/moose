[Mesh]
  [cmg]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${xmin}
    xmax = ${xmax}
    ymin = ${fparse ymin}
    ymax = 0
    nx = 300
    ny = 80
  []
  displacements = 'disp_x disp_y'
[]
