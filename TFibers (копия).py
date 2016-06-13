# coding:utf8
"""
Fibers orientation for left ventrical using 
LDRB algorithm from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3518842/
"""
# Copyright (C) 2015 Zverev Vladimir, Koshelev Anton ( Ural Federal University)
 

from dolfin import *
from math import pi
import numpy as np
from simple_bislerp import bislerp

rad = pi/180
alpha_endo = 40*rad
alpha_epi = 250*rad
beta_endo = 265*rad
beta_epi = 25*rad 

def norm_vector(vector, VectorSpace):
    """
     Function should return normalised vector, 
     but in result will be vector with length more than 1.0
     todo: avoid this situation 
    """
    return project(vector/sqrt(dot(vector, vector)), VectorSpace) 

def axis(grad_psi, grad_phi, ScalarSpace, VectorSpace):
    e_1 = norm_vector(grad_psi, VectorSpace) 
#    plot (e_1)   
    correct_coef = project(dot(e_1, grad_phi), ScalarSpace)
    correct_phi = project(grad_phi - correct_coef*e_1, VectorSpace)
    e_2 = norm_vector(correct_phi, VectorSpace) 
#    plot(e_2)
    e_0 = project(cross(e_1, e_2), VectorSpace) 
#    plot(e_3)
    return (e_0, e_1, e_2)

def orient(Q, alpha, beta, VectorSpace):

    _F, _S, _T = Q
    R_z = as_matrix([[cos(alpha), -sin(alpha), 0],\
                     [sin(alpha),  cos(alpha), 0],\
                     [         0,           0, 1]])

    R_x = as_matrix([[1,          0,         0], 
                     [0,  cos(beta), sin(beta)],
                     [0, -sin(beta), cos(beta)]])
    R_zx = R_z * R_x
    V1 = project(R_zx*_F, VectorSpace)
    V2 = project(R_zx*_S, VectorSpace)
    V3 = project(R_zx*_T, VectorSpace)  
#    plot(V1)  
#    plot(V2)
    plot(V3)

    return (V3, V2, V1)

def bislerp_for_fenics(Qa, Qb, phi, VectorSpace):
    """
    implementation bislerp function which takes into account 
    specific Fenics library

    Original bislerp was described in 
    "Supplementary Material for: A Novel Rule-Based
    Algorithm for Assigning Myocardial Fiber
    Orientation to Computational Heart Models"
    (Bayer, Blake, Plank, Trayanova)

    Input: Qa, Qb - two 3x3 orthogonal matrices
    Input: t ∈ [0, 1] - interpolation factor, 0 → Qa and 1 → Qb
    Input: VectorSpace - space for storing 3d vectors in the Fenics's function
    """
    V1a, V2a, V3a = Qa
    V1b, V2b, V3b = Qb
    
    dof2vertex_scalar = dof_to_vertex_map(phi.function_space())
    phi_values = phi.vector().array()[dof2vertex_scalar]
    vertex2dof_vector = vertex_to_dof_map(VectorSpace)
    dof2vertex = dof_to_vertex_map(VectorSpace)
 
    value_V1a = V1a.vector().array()[dof2vertex]
    value_V2a = V2a.vector().array()[dof2vertex]
    value_V3a = V3a.vector().array()[dof2vertex]

    value_V1b = V1b.vector().array()[dof2vertex]
    value_V2b = V2b.vector().array()[dof2vertex]
    value_V3b = V3b.vector().array()[dof2vertex]

    for n in  np.arange(len(phi_values)):
        s = range(3*n, 3*n + 3)
        Qa_ = np.concatenate((value_V1a[s], value_V2a[s], value_V3a[s])).reshape(3,3)
        Qb_ = np.concatenate((value_V1b[s], value_V2b[s], value_V3b[s])).reshape(3,3)
        R = bislerp(Qa_, Qb_, phi_values[n])
#        print R
        value_V1a[s] = R[:,0]
        value_V2a[s] = R[:,1]  
        value_V3a[s] = R[:,2] 

    V1 = Function(VectorSpace)    
    V1.vector().set_local(value_V1a[vertex2dof_vector])
    V2 = Function(VectorSpace)    
    V2.vector().set_local(value_V2a[vertex2dof_vector])
    V3 = Function(VectorSpace)    
    V3.vector().set_local(value_V3a[vertex2dof_vector])
    return (V1, V2, V3)


def CreateFST(mesh, subdomains, boundaries):
    if mpi_comm_world().rank == 0:
        V = FunctionSpace(mesh, "Lagrange", 1)    
        gdim = mesh.geometry().dim()
        a0 = mesh.coordinates()[:, 0].min()
        a1 = mesh.coordinates()[:, 0].max()
    #    boundaries = FacetFunction('size_t', mesh)
        apex = AutoSubDomain(lambda x: abs(x[0] - a1) < (a1-a0)/200.0)
        apex.mark(boundaries, 7)
        
        """
        legend from lv_from_ply.msh
        ? № Name
        1 1 "ENDORING"
        1 2 "EPIRING"
        2 3 "ENDO"
        2 4 "EPI"
        2 5 "BASE"
        3 6 "MYOCARDIUM"
        """
        bcs_phi = [DirichletBC(V, 1.0, boundaries, 4),
               DirichletBC(V, 0.0, boundaries, 3)]
        bcs_psi = [DirichletBC(V, 1.0, boundaries, 5),
               DirichletBC(V, 0.0, boundaries, 7)]
        
    #    bcs = []
        
        u = TrialFunction(V)
        v = TestFunction(V)
        a = inner(nabla_grad(u), nabla_grad(v)) * dx
        L = Constant(0) * v * dx
        
        phi_epi = Function(V) # for case without RV phi_epi = phi_lv, phi_rv = 0
        solve(a == L, phi_epi, bcs_phi)
        psi = Function(V)
        solve(a == L, psi, bcs_psi)
    #    plot(phi_epi)
    #    plot(psi)
        
        grad_phiepi = grad(phi_epi)
        grad_psi = grad(psi)
    #    plot(grad_phiepi)
    #    plot(grad_psi)

        # phi_rv = 0 hense only the values of the alpha_s and beta_s  at the point 0 is of interest 
        alpha_s = project(alpha_endo, V)
        beta_s = project(beta_endo, V)   

        alpha_w = project(alpha_endo*(1 - phi_epi) + alpha_epi*phi_epi, V)
        beta_w = project(beta_endo*(1 - phi_epi) + beta_epi*phi_epi, V)   

        VectorSpace = VectorFunctionSpace(mesh, "Lagrange", 1)


        # Q_endo = Q_lv in our case, because phi_rv = 0
        Triple_endo = axis(grad_psi, -1*grad_phiepi, V, VectorSpace)
        Q_endo = orient(Triple_endo, alpha_s, beta_s, VectorSpace) 
#        plot(Q_endo[0])
        Triple_epi = axis(grad_psi, grad_phiepi, V, VectorSpace)
        Q_epi = orient(Triple_epi, alpha_w, beta_w, VectorSpace)
#        plot(Q_epi[0])

        (F, S, T) = bislerp_for_fenics (Q_endo, Q_epi, phi_epi, VectorSpace)
#        (F, S, T) = Q_epi
        plot(F)
        interactive()
        F.rename("FiberDir", "tangent vector fiber")
#        File("F_vector.xdmf") << F
#        File("Q_endo_F.xdmf") << Q_endo[0]
#        File("Q_epi_F.xdmf") << Q_epi[0]
#        File("endo_to_epi.xdmf") << phi_epi
        return F

if __name__ == "__main__":
   
    lv_mesh = Mesh("lv_mesh.xml")   
    lv_subdomains = MeshFunction("size_t", lv_mesh, "lv_mesh_physical_region.xml")
    lv_boundaries = MeshFunction("size_t", lv_mesh, "lv_mesh_facet_region.xml")
        
    F =  CreateFST(lv_mesh, lv_subdomains, lv_boundaries)

    hdf = HDF5File(lv_mesh.mpi_comm(), "lv_marked_mesh_with_fibers.h5", "w")
    hdf.write(lv_mesh, "/mesh")
    hdf.write(lv_subdomains, "/subdomains")
    hdf.write(lv_boundaries, "/boundaries")
    hdf.write(F, 'F')
