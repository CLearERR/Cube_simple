# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 19:41:22 2015


"""

import numpy as np

from quaternion_lib.Quaternion import quaternion

def bislerp(Qa, Qb, t):
    """
    implementation bislerp function from
    "Supplementary Material for: A Novel Rule-Based
    Algorithm for Assigning Myocardial Fiber
    Orientation to Computational Heart Models"
    (Bayer, Blake, Plank, Trayanova)

    Input: Qa, Qb - two 3x3 orthogonal matrices
    Input: t ∈ [0, 1] - interpolation factor, 0 → Qa and 1 → Qb
    """
    qA = quaternion(Qa)
    qB = quaternion(Qb)
    unit_quat = [quaternion( 1, 0, 0, 0), 
                 quaternion(-1, 0, 0, 0),
                 quaternion(0,  1, 0, 0), 
                 quaternion(0, -1, 0, 0),
                 quaternion(0, 0,  1, 0), 
                 quaternion(0, 0, -1, 0),
                 quaternion(0, 0, 0,  1), 
                 quaternion(0, 0, 0, -1)]
    qm_list = [qU*qA for qU in unit_quat]
    _test_norm = []
    for qm in qm_list:
        tmp = qm*qB
        _test_norm.append(tmp.norm())
    test_norm = np.array(_test_norm)
    qM = qm_list[np.argmax(test_norm)]
    q = quaternion.slerp_interp(qM, qB, t)
    return q.r()

if __name__ == '__main__':
    from math import pi as PI

#    a = PI
#    ca = np.cos(a)
#    sa = np.sin(a)
#    Ra = np.matrix([[1,  0,    0],
#                    [0,  ca, -sa],
#                    [0,  sa,  ca]]) 
#    print "Ra.shape = {}".format(Ra.shape)  
#    b = PI/2
#    cb = np.cos(b)
#    sb = np.sin(b)
#    Rb = np.matrix([[cb,   0,   sb],
#                    [0,    1,    0],
#                    [-sb,  0,   cb]])  
#    print "Ra:\n", Ra
#    print "Rb:\n", Rb
#    print "bislerp(Ra, Rb, 0):\n", bislerp(Ra, Rb, 0)
#    print "bislerp(Ra, Rb, 1):\n ", bislerp(Ra, Rb, 1)
    np.seterr(all='raise')
    phi = 0.861722369335
    Ra = np.matrix([[-0.5863827,  -0.4420977,  -0.67909161],
                  [-0.14878792,  0.88313441, -0.44664537],
                  [ 0.79677768, -0.16087421, -0.58316854]])

    Rb = np.matrix([[ 0.62120414,  0.44238496,  0.6567119 ],
                    [ 0.6911025,  -0.78690482, -0.03548273],
                    [ 0.50345536,  0.38693638, -0.8514658 ]])
    
    v1 = np.squeeze(np.asarray(Rb[:,0]))
    v2 = np.squeeze(np.asarray(Rb[:,1]))
    v3 = np.squeeze(np.asarray(Rb[:,2]))

    print "check ortonormal Ra", Ra*Ra.T
    print "check ortonormal Rb", Rb*Rb.T
    print "np.dot(v1,v2) = {}".format(np.dot(v1,v2))
    print "np.dot(v2,v3) = {}".format(np.dot(v2,v3)) 
    print "np.dot(v1,v3) = {}".format(np.dot(v1,v3))  

    print "det(Ra) = {}".format(np.linalg.det(Ra))
    print "det(Rb) = {}".format(np.linalg.det(Rb))

    print "bislerp(Ra, Rb, phi):\n", bislerp(Ra, Rb, phi)
