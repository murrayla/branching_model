// Author: Liam Andrew Murray
// email: murrayla@student.unimelb.edu.au
// file: myo_SESM1.geo

// Gmsh project created on Wed Sep 27 12:04:18 2023
SetFactory("OpenCASCADE");

// ==== PARAMETERS ==== //
SARC_R = 0.5;
SARC_L = 2;
MSH_S = 0.1;
// Circle Geom
// Clockwise from Top then Centre
CX = {0, SARC_R, 0, -SARC_R, 0};
CY = {SARC_R, 0, -SARC_R, 0, 0};
// Ellipse Geom
// Clockwise from Top then Centre, Major Axis
// Bottom
EX_B = {0, SARC_R, 0, -SARC_R, 0, 0.1};
EY_B = {0, -SARC_R / 2, -SARC_R, -SARC_R / 2, -SARC_R / 2, -SARC_R / 2};
// Top
EX_T = {0, SARC_R, 0, -SARC_R, 0, 0.1};
EY_T = {SARC_R, SARC_R / 2, 0, SARC_R / 2, SARC_R / 2, SARC_R / 2};
// Z-Disc 0 Params
Z1_N = 1;
ZD0_PN = 5;
Z0 = 0;
// Z-Disc 1 Params
Z2_N = 2;
ZD1_PN = 5;
Z1 = Z0 + 1 * SARC_L;
// Z-Disc 2 Params
Z3_N = 3;
GAP = 0.1;
ZD2_ET_PN = 6;
ZD2_EB_PN = 6;
Z2 = Z0 + 2 * SARC_L;
// Z-Disc 3 Params
Z4_N = 4;
ZD3_ET_PN = 6;
ZD3_EB_PN = 6;
Z3 = Z0 + 3 * SARC_L;
// SPLITS ARE ALONG Y-AXIS, LOOKING DOWN Z-AXIS - to +
POI_N = 0;
ENT_N = 0;
LOO_N = 1;
SUR_N = 1;
VOL_N = 1;

// ==== BEFORE SPLIT [0um] ==== //
// Points and Lines
For i In {1:ZD0_PN}
    Point(i) = {CX[i-1], CY[i-1], Z0, MSH_S};
EndFor
//
For j In {1:(ZD0_PN-1)}
    If (j != (ZD0_PN-1))
        Circle(j) = {j, ZD0_PN, j+1};
    EndIf
    If (j == (ZD0_PN-1))
        Circle(j) = {j, ZD0_PN, (j-ZD0_PN+2)};
    EndIf
EndFor
// Loops and Surfaces
Curve Loop(LOO_N) = {1:(ZD0_PN-1)};
Plane Surface(SUR_N) = {Z1_N};
//
POI_N = POI_N + i - 1;
ENT_N = ENT_N + j - 1;
LOO_N = LOO_N + 2;
SUR_N = SUR_N + 1;

// ==== SPLIT POINT [2um] ==== //
// Points and Lines
For i In {1:ZD1_PN}
    Point(i + POI_N) = {CX[i-1], CY[i-1], Z1, MSH_S};
EndFor
//
For j In {1:ZD1_PN-1}
    If (j != (ZD1_PN-1))
        Circle(j + ENT_N) = {j + i - 1, ZD0_PN + ZD1_PN, j + POI_N + 1};
    EndIf
    If (j == (ZD1_PN-1))
        Circle(j + ENT_N) = {j + i - 1, ZD0_PN + ZD1_PN, j - ZD1_PN + POI_N + 2};
    EndIf
EndFor
//
POI_N = POI_N + i - 1;
ENT_N = ENT_N + j - 1;
//
Line(ENT_N + 1) = {POI_N, POI_N - 3};
Line(ENT_N + 2) = {POI_N, POI_N - 1};
//
ENT_N = ENT_N + 2;
// Loops and Surfaces
Curve Loop(LOO_N) = {ENT_N - 1, ENT_N, ENT_N - 3, ENT_N - 4};
Plane Surface(SUR_N) = {LOO_N};
Curve Loop(LOO_N + 2) = {ENT_N - 1, ENT_N, ENT_N - 2, ENT_N - 5};
Plane Surface(SUR_N + 1) = {LOO_N + 2};
//
LOO_N = LOO_N + 4;
SUR_N = SUR_N + 2;

// ==== CONNECTION [0-2um] ==== //
For k In {1:ZD0_PN - 1}
    BSpline(ENT_N + k) = {k, ZD0_PN + k};
EndFor
//
ENT_N = ENT_N + k - 1;
//
For z In {1:ZD0_PN - 1}
    If (z == 1)
        Curve Loop(LOO_N + 2 * z) = {ZD0_PN - 1, ENT_N - 3, - (ENT_N - 6), - ENT_N};
    EndIf
    If (z != 1)
        Curve Loop(LOO_N + 2 * z) = {ZD0_PN - z, ENT_N - z + 1, - (ENT_N - 5 - z), - (ENT_N - z + 2)};
    EndIf
    Surface(SUR_N + z - 1) = {LOO_N + 2 * z};
EndFor
//
SUR_N = SUR_N + z - 1; 
LOO_N = LOO_N + ZD0_PN * 2; 
//
Surface Loop(VOL_N) = {1:(SUR_N - 1)};
Volume(VOL_N) = {VOL_N};
//
VOL_N = VOL_N + 1;

// ==== SPLIT END [4um] ==== //
// Points and Lines
// Bottom Ellipse
For i In {1:ZD2_EB_PN}
    Point(i + POI_N) = {EX_B[i-1], EY_B[i-1] - GAP, Z2, MSH_S}; 
EndFor
//
For j In {1:ZD2_EB_PN - 2}
    If (j == (ZD1_PN-1))
        Ellipse(j + ENT_N) = {POI_N + 4, POI_N + 5, POI_N + 6, POI_N + 1};
    Else
        Ellipse(j + ENT_N) = {POI_N + j, POI_N + 5, POI_N + 6, POI_N + j + 1};
    EndIf
EndFor
//
Curve Loop(LOO_N) = {ENT_N + 1, ENT_N + 2, ENT_N + 3, ENT_N + 4};
Surface(SUR_N) = {LOO_N};
//
POI_N = POI_N + i - 1;
ENT_N = ENT_N + j - 1;
LOO_N = LOO_N + 2;
SUR_N = SUR_N + 1;
// Top Ellipse
For i In {1:ZD2_ET_PN}
    Point(i + POI_N) = {EX_T[i-1], EY_T[i-1] + GAP, Z2, MSH_S}; 
EndFor
//
For j In {1:ZD2_ET_PN - 2}
    If (j == (ZD1_PN-1))
        Ellipse(j + ENT_N) = {POI_N + 4, POI_N + 5, POI_N + 6, POI_N + 1};
    Else
        Ellipse(j + ENT_N) = {POI_N + j, POI_N + 5, POI_N + 6, POI_N + j + 1};
    EndIf
EndFor
//
Curve Loop(LOO_N) = {ENT_N + 1, ENT_N + 2, ENT_N + 3, ENT_N + 4};
Surface(SUR_N) = {LOO_N};
//
POI_N = POI_N + i - 1;
ENT_N = ENT_N + j - 1;
LOO_N = LOO_N + 1;
SUR_N = SUR_N + 1;

// ==== CONNECTION [2-4um] ==== //
// Bottom Ellipse
For k In {1:ZD2_EB_PN - 2}
    If (k == 1)
        BSpline(ENT_N + k) = {POI_N - 13 + k, POI_N - 12 + k};
    Else
        BSpline(ENT_N + k) = {POI_N - 17 + k, POI_N - 12 + k};
    EndIf
EndFor
//
ENT_N = ENT_N + k - 1;
//
For z In {1:ZD2_EB_PN - 2}
    If (z == 1)
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 15, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z - 2)};
    ElseIf (z == 4)
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 13, ENT_N - z + 2, - (ENT_N - z - 7), - (ENT_N - z + 1)};
    Else
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 17, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z + 2)};
    EndIf
    Surface(SUR_N + z - 1) = {LOO_N + 2 * z};
EndFor
//
SUR_N = SUR_N + z - 1; 
LOO_N = LOO_N + ZD0_PN * 2; 
//
Surface Loop(VOL_N) = {SUR_N - 12, SUR_N - 6, (SUR_N - 4):(SUR_N - 1)};
Volume(VOL_N) = {VOL_N};
//
VOL_N = VOL_N + 1;
// Top Ellipse
For k In {1:ZD2_ET_PN - 2}
    If (k >= 3)
        BSpline(ENT_N + k) = {POI_N - 9 - k, POI_N - 6 + k};
    Else
        BSpline(ENT_N + k) = {POI_N - 17 + k, POI_N - 6 + k};
    EndIf
EndFor
//
ENT_N = ENT_N + k - 1;
//
For z In {1:ZD2_ET_PN - 2}
    If (z == 1)
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 21, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z - 2)};
    ElseIf (z == 2)
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 18, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z + 2)};
    ElseIf (z == 3)
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 18, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z + 2)};
    Else
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 21, ENT_N - z + 2, - (ENT_N - z - 7), - (ENT_N - z + 1)};
    EndIf
    Surface(SUR_N + z - 1) = {LOO_N + 2 * z};
EndFor
//
SUR_N = SUR_N + z - 1; 
LOO_N = LOO_N + ZD0_PN * 2; 
//
Surface Loop(VOL_N) = {SUR_N - 15, SUR_N - 9, (SUR_N - 4):(SUR_N - 1)};
Volume(VOL_N) = {VOL_N};
//
VOL_N = VOL_N + 1;

// ==== END [6um] ==== //
// Points and Lines
// Bottom Ellipse
For i In {1:ZD3_EB_PN}
    Point(i + POI_N) = {EX_B[i-1], EY_B[i-1] - GAP, Z3, MSH_S}; 
EndFor
//
For j In {1:ZD3_EB_PN - 2}
    If (j == (ZD1_PN-1))
        Ellipse(j + ENT_N) = {POI_N + 4, POI_N + 5, POI_N + 6, POI_N + 1};
    Else
        Ellipse(j + ENT_N) = {POI_N + j, POI_N + 5, POI_N + 6, POI_N + j + 1};
    EndIf
EndFor
//
Curve Loop(LOO_N) = {ENT_N + 1, ENT_N + 2, ENT_N + 3, ENT_N + 4};
Surface(SUR_N) = {LOO_N};
//
POI_N = POI_N + i - 1;
ENT_N = ENT_N + j - 1;
LOO_N = LOO_N + 2;
SUR_N = SUR_N + 1;
// Top Ellipse
For i In {1:ZD3_ET_PN}
    Point(i + POI_N) = {EX_T[i-1], EY_T[i-1] + GAP, Z3, MSH_S}; 
EndFor
//
For j In {1:ZD2_ET_PN - 2}
    If (j == (ZD1_PN-1))
        Ellipse(j + ENT_N) = {POI_N + 4, POI_N + 5, POI_N + 6, POI_N + 1};
    Else
        Ellipse(j + ENT_N) = {POI_N + j, POI_N + 5, POI_N + 6, POI_N + j + 1};
    EndIf
EndFor
//
Curve Loop(LOO_N) = {ENT_N + 1, ENT_N + 2, ENT_N + 3, ENT_N + 4};
Surface(SUR_N) = {LOO_N};
//
POI_N = POI_N + i - 1;
ENT_N = ENT_N + j - 1;
LOO_N = LOO_N + 1;
SUR_N = SUR_N + 1;

// ==== CONNECTION [4-6um] ==== //
// Bottom Ellipse
For k In {1:ZD2_EB_PN - 2}
        BSpline(ENT_N + k) = {POI_N - 24 + k, POI_N - 12 + k};
EndFor
//
ENT_N = ENT_N + k - 1;
//
For z In {1:ZD2_EB_PN - 2}
    If (z == 1)
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 23, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z - 2)};
    Else
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 23, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z + 2)};
    EndIf
    Surface(SUR_N + z - 1) = {LOO_N + 2 * z};
EndFor
//
SUR_N = SUR_N + z - 1; 
LOO_N = LOO_N + ZD0_PN * 2; 
//
Surface Loop(VOL_N) = {SUR_N - 16, SUR_N - 6, (SUR_N - 4):(SUR_N - 1)};
Volume(VOL_N) = {VOL_N};
//
VOL_N = VOL_N + 1;
// Top Ellipse
For k In {1:ZD2_ET_PN - 2}
    BSpline(ENT_N + k) = {POI_N - 18 + k, POI_N - 6 + k};
EndFor
//
ENT_N = ENT_N + k - 1;
//
For z In {1:ZD2_ET_PN - 2}
    If (z == 1)
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 23, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z - 2)};
    Else
        Curve Loop(LOO_N + 2 * z) = {ENT_N - z - 23, ENT_N - z + 1, - (ENT_N - z - 7), - (ENT_N - z + 2)};
    EndIf
    Surface(SUR_N + z - 1) = {LOO_N + 2 * z};
EndFor
//
SUR_N = SUR_N + z - 1; 
LOO_N = LOO_N + ZD0_PN * 2; 
//
Surface Loop(VOL_N) = {SUR_N - 19, SUR_N - 9, (SUR_N - 4):(SUR_N - 1)};
Volume(VOL_N) = {VOL_N};
//
VOL_N = VOL_N + 1;
