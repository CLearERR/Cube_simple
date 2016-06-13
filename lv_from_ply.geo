// meshing options
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.Lloyd = 1;
Mesh.CharacteristicLengthMin = 3.0;
Mesh.CharacteristicLengthMax = 4.0;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.RemeshParametrization = 7;
Mesh.SurfaceFaces = 1;
Mesh.CharacteristicLengthFactor = 0.5;
Mesh.RemeshAlgorithm = 1;

// load the surfaces
Merge "Endo.ply";
Merge "Epi.ply";

CreateTopology;

ll[] = Line "*";
L_LV_base = newl; Compound Line(L_LV_base) = ll[1];
L_epi_base = newl; Compound Line(L_epi_base) = ll[0];
Physical Line("ENDORING") = { L_LV_base };
Physical Line("EPIRING") = { L_epi_base };

ss[] = Surface "*";
S_LV = news; Compound Surface(S_LV) = ss[0];
S_epi = news; Compound Surface(S_epi) = ss[1];
Physical Surface("ENDO") = { S_LV };
Physical Surface("EPI") = { S_epi };

LL_base = newll; 
Line Loop(LL_base) = { L_LV_base, L_epi_base };
S_base = news; Plane Surface(S_base) = { LL_base };
Physical Surface("BASE") = { S_base };

SL_wall = newsl; 
Surface Loop(SL_wall) = { S_LV, S_epi, S_base };
V_wall = newv; Volume(V_wall) = { SL_wall };
Physical Volume("MYOCARDIUM") = { V_wall };

Physical Surface(12) = {11};
