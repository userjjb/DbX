from dolfin import *

# the L-shape mesh can be find at
# http://fenicsproject.org/pub/data/meshes/
mesh = Mesh("lshape.xml.gz")
mesh.coordinates()[:] = (mesh.coordinates()[:] -0.5)*2.0

def D_boundary(x, on_boundary):
    return on_boundary

#V = FunctionSpace(mesh,'CG',1)
W = VectorElement ("CG",mesh.ufl_cell(), 2 )
Q = FiniteElement ("CG",mesh.ufl_cell(), 1 )
V = FunctionSpace(mesh, W * Q)
uh = TrialFunction(V)
vh = TestFunction(V)

ue = Expression("pow(x[0]*x[0]+x[1]*x[1],1.0/3.0)", element=V.ufl_element())
f = Expression("-4.0/9.0*pow(x[0]*x[0]+x[1]*x[1],-2.0/3.0)", element=V.ufl_element())

a = inner(grad(uh),grad(vh))*dx
L = f*vh*dx

A = assemble(a)
b = assemble(L)
bc = DirichletBC(V, ue, D_boundary)
bc.apply(A,b)

#---------------------------------------------------
#solver = KrylovSolver("gmres", "ml_amg")
#solver.parameters["absolute_tolerance"] = 1E-13
#solver.parameters["relative_tolerance"] = 1E-13
#solver.parameters["maximum_iterations"] = 1000000
#solver.parameters["nonzero_initial_guess"] = True
#--------------------------------------------
uh = Function(V)
#uh.vector()[:] = uh.vector()[:]*0.8
solve(A, uh.vector(), b)

EL2 = errornorm(ue,uh,"L2", degree_rise = 1)
EH1 = errornorm(ue,uh,"H1", degree_rise = 1)

print (EL2, EH1)