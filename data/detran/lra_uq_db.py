from detran import *

def get_input() :
  inp = utilities.InputDB.Create()
  inp.put_int("number_groups",                  2)
  inp.put_int("dimension",                      2)
  inp.put_str("equation",                       "diffusion")
  inp.put_str("bc_west",                        "reflect")
  inp.put_str("bc_east",                        "vacuum")
  inp.put_str("bc_south",                       "reflect")
  inp.put_str("bc_north",                       "vacuum")
  inp.put_int("bc_zero_flux",                   0)
  inp.put_int("quad_number_polar_octant",       3)
  inp.put_int("quad_number_azimuth_octant",     3)
  inp.put_str("eigen_solver",                   "arnoldi")
  inp.put_dbl("eigen_tolerance",                1e-12)
  inp.put_int("eigen_max_iters",                1000)
  inp.put_str("outer_solver",                   "GMRES")
  inp.put_dbl("outer_tolerance",                1e-12)
  inp.put_int("outer_max_iters",                1000)
  inp.put_int("outer_print_level",              0)
  inp.put_int("outer_krylov_group_cutoff",      0)
  inp.put_str("outer_pc_type",                  "mgdsa")
  inp.put_str("inner_solver",                   "GMRES")
  inp.put_dbl("inner_tolerance",                1e-12)
  inp.put_int("inner_max_iters",                1000)
  inp.put_int("inner_print_level",              0)
  inp.put_str("inner_pc_type",                  "DSA")
  # gmres parameters
  db = utilities.InputDB.Create("callow_db")
  #db.put_dbl("linear_solver_atol",              0.0);
  db.put_dbl("linear_solver_rtol",              1e-12);
  db.put_str("linear_solver_type",              "gmres");
  db.put_int("linear_solver_maxit",             1000);
  db.put_int("linear_solver_gmres_restart",     30);
  db.put_int("linear_solver_monitor_level",     0);
  db.put_str("pc_type",                         "jacobi");
  db.put_str("petsc_pc_type",                   "lu");
  db.put_int("petsc_pc_factor_levels",          3);
  db.put_str("eigen_solver_type",               "gd");
  db.put_int("eigen_solver_maxit",              1000);
  db.put_int("eigen_solver_monitor_level",      1);
  db.put_dbl("eigen_solver_tol",                1.0e-9)
  inp.put_spdb("inner_solver_db", db)
  inp.put_spdb("inner_pc_db", db)
  inp.put_spdb("outer_solver_db", db) 
  inp.put_spdb("eigen_solver_db", db)
  inp.put_int("ts_max_steps",                   10000)
  inp.put_int("ts_scheme",                      Time2D.BDF2)
  inp.put_int("ts_output",                      0)
  inp.put_dbl("ts_step_size",                   0.01)
  inp.put_dbl("ts_final_time",                  3.0)
  #inp.put_int("ts_no_extrapolation",            1)
  inp.put_int("ts_max_iters",                   10)
  inp.put_dbl("ts_tolerance",                   1.0e-5)
  #
  preconditioner_db = utilities.InputDB.Create("preconditioner_db")
  preconditioner_db.put_dbl("linear_solver_atol",              0.0);
  preconditioner_db.put_dbl("linear_solver_rtol",              0.1);
  preconditioner_db.put_str("linear_solver_type",              "gmres");
  preconditioner_db.put_int("linear_solver_maxit",             5000);
  preconditioner_db.put_int("linear_solver_gmres_restart",     30);
  preconditioner_db.put_int("linear_solver_monitor_level",     0);
  preconditioner_db.put_str("pc_type",                         "ilu0");
  preconditioner_db.put_str("petsc_pc_type",                   "ilu");
  preconditioner_db.put_int("petsc_pc_factor_levels",          2);
  #
  return inp

