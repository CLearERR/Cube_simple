# coding:utf8

from dolfin import *
import numpy as np
from simple_bislerp import bislerp

def test2():
    mesh = BoxMesh(0,1,0,1,0,1, 8,4,4)
    print "num_vertices() = {}".format(mesh.num_vertices())
    ScalarSpace = FunctionSpace(mesh, "Lagrange", 1)
    S = interpolate(Expression("x[0]"), ScalarSpace)
    dof2vertex_scalar = dof_to_vertex_map(ScalarSpace)
    S_values = S.vector().array()[dof2vertex_scalar]


    VectorSpace = VectorFunctionSpace(mesh, "Lagrange", 1)
    
    e1 = Constant((1.0, 0.0, 0.0))
    e2 = Constant((0.0, 1.0, 0.0))
    e3 = Constant((0.0, 0.0, 1.0))    
    Q1 = interpolate(e1, VectorSpace)
    Q2 = interpolate(e2, VectorSpace)
    Q3 = interpolate(e3, VectorSpace)

    e1_s = Constant((-1.0, 0.0, 0.0))
    e2_s = Constant((0.0, -1.0, 0.0))
    e3_s = Constant((0.0, 0.0, 1.0))    
    Q1_s = interpolate(e1_s, VectorSpace)
    Q2_s = interpolate(e2_s, VectorSpace)
    Q3_s = interpolate(e3_s, VectorSpace)
#    plot(Q1)
#    plot(Q2)
#    plot(Q3)
    interactive()
    vertex2dof = vertex_to_dof_map(VectorSpace)
    dof2vertex = dof_to_vertex_map(VectorSpace)

    
    value_q1 = Q1.vector().array()[dof2vertex]
    value_q2 = Q2.vector().array()[dof2vertex]
    value_q3 = Q3.vector().array()[dof2vertex]

    value_q1_s = Q1_s.vector().array()[dof2vertex]
    value_q2_s = Q2_s.vector().array()[dof2vertex]
    value_q3_s = Q3_s.vector().array()[dof2vertex]

    array_dim = len(value_q1)
    print "len(value_q1) = {}".format(array_dim)
    
    n = 0

    # we use vector from 3d space therefore we multiply index of node by 3
    # n = 0 it is not first node, all rather comlicate dur internal dof numbering
#    value_q1[3*n] = 2
#    value_q1[3*n+1] = 2
#    value_q1[3*n+2] = 2
    
#    value_q2[3*n] = 3
#    value_q2[3*n+1] = 3
#    value_q2[3*n+2] = 3
#    
#    value_q3[3*n] = 4
#    value_q3[3*n+1] = 4
#    value_q3[3*n+2] = 4

#    print value_q2    

#    tmp1 = np.concatenate((value_q1, value_q2, value_q3))
#    tmp2 = tmp1.reshape(3,len(value_q1))
#    joint_matrix = tmp2.T
    
    for n in  np.arange(len(S_values)):
        s = range(3*n, 3*n + 3)
        Qa = np.concatenate((value_q1[s], value_q2[s], value_q3[s])).reshape(3,3)
        Qb = np.concatenate((value_q1_s[s], value_q2_s[s], value_q3_s[s])).reshape(3,3)
        R = bislerp(Qa, Qb, S_values[n])
#        print R
        value_q1[s] = R[:,0]
        value_q2[s] = R[:,1]  
        value_q3[s] = R[:,2] 
        
    Q1.vector().set_local(value_q1[vertex2dof])
    Q1file = File("Q1test.pvd")
    Q1file << Q1
    Sfile = File("S_forQ1test.pvd")
    Sfile << S
#    plot(Q1)
#    interactive()
   

def test1():
    M=UnitSquareMesh(2,2)

    TS=TensorFunctionSpace(M, "CG", 1)

    vd=vertex_to_dof_map(TS)
    dv=dof_to_vertex_map(TS)

    F=Function(TS)

    F=interpolate(Expression((('2*x[0]','x[1]'),('1.0','2.0'))),TS)

    ve=F.vector().array()[dv]

    #Change the values for node 2, for instance
    n=2

    ve[n*4]=0
    ve[n*4+1]=1
    ve[n*4+2]=2
    ve[n*4+3]=3

    F.vector().set_local(ve[vd])

    #plot one component to see what happened
    plot(F[0,0])
    interactive()
    
if __name__ == "__main__":
    test2()
