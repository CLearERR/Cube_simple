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
    vector  = Function(VectorSpace)
    Function should return normalised vector, 
    but in result will be vector with length more than 1.0
    todo: avoid this situation 
    
    u1 = Function(V) 
    u1.vector()[:] = u0.vector() + s*delta_u.vector()
    """
    

    vertex2dof_vector = vertex_to_dof_map(VectorSpace)
    dof2vertex_vector = dof_to_vertex_map(VectorSpace)
 
    values = vector.vector().array()[dof2vertex_vector]
    for n in  np.arange(len(values)/3):
        s = range(3*n, 3*n + 3)
        v3d = values[s]
        values[s] = v3d/np.linalg.norm(v3d)
#        print np.linalg.norm(v3d), np.linalg.norm(values[s])
        
    tmp = Function(VectorSpace)    
    tmp.vector().set_local(values[vertex2dof_vector])
    
    #return project(vector/norm_, VectorSpace)
    return tmp
        

def axis(grad_psi, grad_phi, ScalarSpace, VectorSpace):
    """ 
    e0 represents circumference direction
    e1 for apicobasal direction
    e2 for transmural direction
    """
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
    """
    ei = (ei0, ei1, ei2)^T, where i = 1,2,3
        e00 | e10 | e20             e00 | e01 | e02    
    Q = e01 | e11 | e21;       Q' = e10 | e11 | e12 = Q^T;   
        e02 | e12 | e22             e20 | e21 | e22   
 
           cos(a) | -sin(a) | 0               1 |      0 |      0
  R_e2_a = sin(a) | cos(a)  | 0 ; R_F_b =     0 |  cos(b)| sin(b)           
                0 |     0   | 1               0 | -sin(b)| cos(b)

    """
    e0,e1,e2 = Q
    Q1 =  as_matrix([[e0[0], e1[0], e2[0]],
                     [e0[1], e1[1], e2[1]],
                     [e0[2], e1[2], e2[2]]])                             
    ca, sa = cos(alpha), sin(alpha)
    R_e2_a = as_matrix([[ca, -sa, 0],
                        [sa,  ca, 0],
                        [ 0,   0, 1]])                             
                        
    Q2 = Q1*R_e2_a*Q1.T
    cb, sb = cos(beta), sin(beta)
    R_F_b = as_matrix([[0, -sb, cb],
                       [0,  cb, sb],           
                       [1,   0,  0]])
                       
    M = Q2*R_F_b*Q2.T   
    
#    E1 = as_vector([M[0,0],M[0,1]], M[0,2])
#    E2 = as_vector([M[1,0],M[1,1]], M[1,2])
#    E3 = as_vector([M[2,0],M[2,1]], M[2,2])

    E1 = as_vector([1, 0, 0])
    E2 = as_vector([0, 1, 0])
    E3 = as_vector([0, 0, 1])                

    V1 = project(M*E1, VectorSpace)
    V2 = project(M*E2, VectorSpace)
    V3 = project(M*E3, VectorSpace)  
   
#    _F, _S, _T = Q

#    F_new = as_vector([cos(alpha)*_F[0] + sin(alpha)*_S[0],\
#                       cos(alpha)*_F[1] + sin(alpha)*_S[1],\
#                       cos(alpha)*_F[2] + sin(alpha)*_S[2]])
#    
#    S_new = as_vector([cos(beta)*(cos(alpha)*_S[0] - sin(alpha)*_F[0]) - sin(beta)*_T[0],\
#                       cos(beta)*(cos(alpha)*_S[1] - sin(alpha)*_F[1]) - sin(beta)*_T[1],\
#                       cos(beta)*(cos(alpha)*_S[2] - sin(alpha)*_F[2]) - sin(beta)*_T[2] ])

#    T_new = as_vector([sin(beta)*(cos(alpha)*_S[0] - sin(alpha)*_F[0]) + cos(beta)*_T[0],\
#                       sin(beta)*(cos(alpha)*_S[1] - sin(alpha)*_F[1]) + cos(beta)*_T[1],\
#                       sin(beta)*(cos(alpha)*_S[2] - sin(alpha)*_F[2]) + cos(beta)*_T[2] ]) 

#    V1 = project(F_new, VectorSpace)
#    V2 = project(S_new, VectorSpace)
#    V3 = project(T_new, VectorSpace)  
#    plot(V1)  
#    plot(V2)
#    plot(V3)

    return (V1, V2, V3)

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

    np.seterr(all='raise')
    for n in  np.arange(len(phi_values))[0:4]:
        s = range(3*n, 3*n + 3)
        Qa_ = np.concatenate((value_V1a[s], value_V2a[s], value_V3a[s])).reshape(3,3)
        Qb_ = np.concatenate((value_V1b[s], value_V2b[s], value_V3b[s])).reshape(3,3)
#        R = bislerp(Qa_, Qb_, phi_values[n])
        try:
            R = bislerp(Qa_, Qb_, phi_values[n])
        except FloatingPointError:
            print  "n = {}, phi = {}\n".format(n, phi_values[n]),  
            print Qa_ 
            print Qb_ 
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
        a0 = mesh.coordinates()[:, 2].min()
        a1 = mesh.coordinates()[:, 2].max()
    #    boundaries = FacetFunction('size_t', mesh)
        apex = AutoSubDomain(lambda x: x[2]  < (a1-a0)/50)
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
        bcs_phi_epi = [DirichletBC(V, 1.0, boundaries, 4),
                       DirichletBC(V, 0.0, boundaries, 3)]
        bcs_phi_lv = [DirichletBC(V, 1.0, boundaries, 3),
                      DirichletBC(V, 0.0, boundaries, 4)]
        bcs_psi = [DirichletBC(V, 1.0, boundaries, 5),
                   DirichletBC(V, 0.0, boundaries, 7)]
                
        u = TrialFunction(V)
        v = TestFunction(V)
        a = inner(nabla_grad(u), nabla_grad(v)) * dx
        L = Constant(0) * v * dx
        
        # for case without RV dSigma_lv = dSigma_endo, dSigma_rv =0,  phi_rv = 0 
         
        phi_epi = Function(V)
        solve(a == L, phi_epi, bcs_phi_epi )
        phi_lv = Function(V)
        solve(a == L, phi_lv, bcs_phi_lv)
        psi = Function(V)
        solve(a == L, psi, bcs_psi)
#        plot(phi_epi)
#        plot(psi)
        
        VectorSpace = VectorFunctionSpace(mesh, "Lagrange", 1)
        grad_phiepi = project(grad(phi_epi), VectorSpace)
        negative_grad_philv = project( grad(-1*phi_lv), VectorSpace)
        grad_psi = project(grad(psi), VectorSpace) 

#        plot(negative_grad_philv)

        # phi_rv = 0 hense only the values of the alpha_s and beta_s  at the point 0 is of interest 
        alpha_s = project(alpha_endo, V)
        beta_s = project(beta_endo, V)   

        alpha_w = project(alpha_endo*(1 - phi_epi) + alpha_epi*phi_epi, V)
        beta_w = project(beta_endo*(1 - phi_epi) + beta_epi*phi_epi, V)

        # Q_endo = Q_lv in our case, because phi_rv = 0

        Triple_endo = axis(grad_psi, negative_grad_philv, V, VectorSpace)
#        plot(Triple_endo[0])
        Q_endo = orient(Triple_endo, alpha_s, beta_s, VectorSpace) 
        plot(Q_endo[0])

        Triple_epi = axis(grad_psi, grad_phiepi, V, VectorSpace)
#        plot(Triple_epi[0])
        Q_epi = orient(Triple_epi, alpha_w, beta_w, VectorSpace)
        plot(Q_epi[0])

        (F, S, T) = bislerp_for_fenics (Q_endo, Q_epi, phi_epi, VectorSpace)
##        (F, S, T) = Q_epi
        plot(F)
        interactive()
#        F.rename("FiberDir", "tangent vector fiber")
#        File("F_vector.xdmf") << F
#        File("Q_endo_F.xdmf") << Q_endo[0]
#        File("Q_epi_F.xdmf") << Q_epi[0]
#        File("endo_to_epi.xdmf") << phi_epi
#        return F

if __name__ == "__main__":
    parameters['form_compiler']['representation'] = 'quadrature'

    lv_mesh = Mesh("lv_mesh_sym.xml")   
    lv_subdomains = MeshFunction("size_t", lv_mesh, "lv_mesh_physical_region_sym.xml")
    lv_boundaries = MeshFunction("size_t", lv_mesh, "lv_mesh_facet_region_sym.xml")
    
    parameters["form_compiler"]["quadrature_degree"] = 4
    F = CreateFST(lv_mesh, lv_subdomains, lv_boundaries)

#    hdf = HDF5File(lv_mesh.mpi_comm(), "lv_marked_mesh_with_fibers.h5", "w")
#    hdf.write(lv_mesh, "/mesh")
#    hdf.write(lv_subdomains, "/subdomains")
#    hdf.write(lv_boundaries, "/boundaries")
#    hdf.write(F, 'F')
