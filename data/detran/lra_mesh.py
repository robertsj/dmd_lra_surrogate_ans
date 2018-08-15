from detran import *

def get_mesh_2D(fmm, nodal=False) :
  """ Return the 2-D LRA Mesh
  """

  # uniformly-spaced coarse meshes
  cm = np.linspace(0.0, 165.0, 12)
  # fine mesh per coarse mesh (1x1 = 15x15 cm)
  fm = vec_int(11, fmm)

  # unique coarse mesh map
  cmm = [ 1 , 0 , 0 , 0 , 0 , 1 , 1 , 2 , 2 , 4 , 4,
          0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 4 , 4,
          0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 4 , 4,
          0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 4 , 4,
          0 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 4 , 4,
          1 , 0 , 0 , 0 , 0 , 1 , 1 , 5 , 5 , 4 , 4,
          1 , 0 , 0 , 0 , 0 , 1 , 1 , 5 , 5 , 4 , 4,
          2 , 2 , 2 , 2 , 2 , 2 , 2 , 3 , 4 , 4 , 4,
          2 , 2 , 2 , 2 , 2 , 2 , 2 , 4 , 4 , 4 , 4,
          4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4,
          4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4]
  # bundle map [9,9,9,9,9,9,9,8,7]
  bm =  [ 0,  1,  2,  3,  4,  5,  6,  7,  8, 78, 78,
          9, 10, 11, 12, 13, 14, 15, 16, 17, 78, 78,
         18, 19, 20, 21, 22, 23, 24, 25, 26, 78, 78,
         27, 28, 29, 30, 31, 32, 33, 34, 35, 78, 78,
         36, 37, 38, 39, 40, 41, 42, 43, 44, 78, 78,
         45, 46, 47, 48, 49, 50, 51, 52, 53, 78, 78,
         54, 55, 56, 57, 58, 59, 60, 61, 62, 78, 78,
         63, 64, 65, 66, 67, 68, 69, 70, 78, 78, 78,
         71, 72, 73, 74, 75, 76, 77, 78, 78, 78, 78,
         78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
         78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78]
  
  # material map is defaulted to the cmm

  mesh = Mesh2D.Create(fm, fm, cm, cm, cmm)
  mesh.add_coarse_mesh_map("ASSEMBLY",    bm)
  mesh.add_coarse_mesh_map("SUBASSEMBLY", bm)
  mesh.add_coarse_mesh_map("COARSEMESH",  cmm)

  if not nodal :
    # fine mesh material map
    mm = vec_int(mesh.number_cells(), 0)
    for i in range(0, len(mm)) :
      mm[i]  = i
    mesh.add_mesh_map("MATERIAL", mm)

  return mesh, cmm


