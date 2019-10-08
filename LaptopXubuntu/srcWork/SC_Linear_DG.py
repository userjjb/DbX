import numpy as np

#Solve scalar conservation law of form:
#\frac{\partial u, \partial t} + \frac{\partial f(u), \partial x} = 0
#where u is the solution, and f(u) is the "flux" of u

#We solve with DG by first expressing the PDE in a weak form:
#\int \frac{\partial u, \partial t} \phi + \int \frac{\partial f(u), \partial x} \phi = 0
#where \phi is a test function

#We approximate the true solution (which lives in an infinite vector space)
#with a finite space, as a basis choose the set of polynomials of order N or
#less; so u is:
#u = \sum_{i=0}^N a_i \psi_i(x)
#where \psi_i is the ith order polynomial base, and a_i is the ith basis weight

#If we chose our test function space to be the same as out solution space, and
#choose to use the same basis (i.e. \phi_i = \psi_j if i=j) then we have now
#specified all the discrete approximations to the original weak form (note that
#we must satisfy the weak form for *each* of the test functions /phi_j)

#One last detail: to reduce smoothness requirements on f(u) integrate the second
#integral by parts:
#\int_a^b \frac{\partial f(u), \partial x} \phi = 
#   (u* \phi) \rvert_a^b - \int_a^b f(u) \frac{\partial \phi, \partial x}
#to evaluate u* on the element boundaries (where it is multiply defined) we use
#a numerical flux function to determine the value of u*. Let's choose the upwind
#flux function where u*(x_p) is the u(x_p) of the element that is "upwind" of the
#point x_p