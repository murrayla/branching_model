// Author: Liam Andrew Murray
// email: murrayla@student.unimelb.edu.au
// file: myo_SESM1.geo

// Gmsh project created on Wed Sep 27 12:04:18 2023
SetFactory("OpenCASCADE");

// ==== BEFORE SPLIT [0um] ==== //
// Points and Lines
Point(1) = {0, 0, 0, 1};
Point(2) = {0.5, 0, 0, 1};
Point(4) = {0, 0.5, 0, 1};
Point(5) = {-0.5, 0, 0, 1};
Point(6) = {0, -0.5, 0, 1};
Circle(1) = {5, 1, 6};
Circle(2) = {6, 1, 2};
Circle(3) = {2, 1, 4};
Circle(4) = {4, 1, 5};
// Loops and Surfaces
Curve Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};

// ==== SPLIT POINT [2um] ==== //
// Points and Lines
Point(7) = {0, 0, 2, 1};
Point(8) = {0.5, 0, 2, 1};
Point(9) = {0, 0.5, 2, 1};
Point(10) = {-0.5, 0, 2, 1};
Point(11) = {0, -0.5, 2, 1};
Circle(5) = {10, 7, 11};
Circle(6) = {11, 7, 8};
Circle(7) = {8, 7, 9};
Circle(8) = {9, 7, 10};
// Loops and Surfaces
Curve Loop(3) = {5, 6, 7, 8};
Surface(2) = {3};

// ==== CONNECTION [0-2um] ==== //
BSpline(9) = {5, 10};
BSpline(10) = {6, 11};
BSpline(11) = {2, 8};
BSpline(12) = {4, 9};
Curve Loop(5) = {1, 10, -5, -9};
Surface(3) = {5};
Curve Loop(7) = {2, 11, -6, -10};
Surface(4) = {7};
Curve Loop(9) = {7, -12, -3, 11};
Surface(5) = {9};
Curve Loop(11) = {8, -9, -4, 12};
Surface(6) = {11};
Surface Loop(1) = {3, 6, 5, 4, 2, 1};
Volume(1) = {1};

// ==== SPLIT END [4um] ==== //
// Points and Lines
Point(12) = {0, 0, 4, 1.0};
// Ellipse 1
Point(13) = {0, -0.5, 4, 1.0};
Point(14) = {-0.5, -0.25, 4, 1.0};
Point(15) = {0.5, -0.25, 4, 1.0};
Point(16) = {0, -0.25, 4, 1.0};
Line(28) = {7, 10};
Line(29) = {7, 8};
// Ellipse 2
Point(17) = {0, 0.5, 4, 1.0};
Point(18) = {-0.5, 0.25, 4, 1.0};
Point(19) = {0.5, 0.25, 4, 1.0};
Point(20) = {0, 0.25, 4, 1.0};
// Major Axis Points
Point(21) = {0.1, -0.25, 4, 1.0};
Point(22) = {0.1, 0.25, 4, 1.0};
// Ellipse 1
Ellipse(13) = {14, 16, 21, 13};
Ellipse(14) = {13, 16, 21, 15};
Ellipse(15) = {15, 16, 21, 12};
Ellipse(16) = {12, 16, 21, 14};
// Ellipse 2
Ellipse(17) = {18, 20, 22, 12};
Ellipse(18) = {12, 20, 22, 19};
Ellipse(19) = {19, 20, 22, 17};
Ellipse(20) = {17, 20, 22, 18};
// Loops and Surfaces
Curve Loop(29) = {29, 7, 8, -28};
Surface(15) = {29};
Curve Loop(31) = {28, 5, 6, -29};
Surface(16) = {31};
Curve Loop(41) = {20, 17, 18, 19};
Surface(17) = {41};
Curve Loop(43) = {16, 13, 14, 15};
Surface(18) = {43};

// ==== CONNECTION [2-4um] ==== //
BSpline(21) = {7, 12};
BSpline(22) = {10, 14};
BSpline(23) = {11, 13};
BSpline(24) = {8, 15};
BSpline(25) = {10, 18};
BSpline(26) = {9, 17};
BSpline(27) = {8, 19};
Curve Loop(13) = {22, 13, -23, -5};
Surface(7) = {13};
Curve Loop(15) = {23, 14, -24, -6};
Surface(8) = {15};
Curve Loop(17) = {15, -21, 29, 24};
Surface(9) = {17};
Curve Loop(19) = {21, 16, -22, -28};
Surface(10) = {19};
Curve Loop(21) = {8, 25, -20, -26};
Surface(11) = {21};
Curve Loop(23) = {7, 26, -19, -27};
Surface(12) = {23};
Curve Loop(25) = {18, -27, -29, 21};
Surface(13) = {25};
Curve Loop(27) = {17, -21, 28, 25};
Surface(14) = {27};
Surface Loop(2) = {12, 11, 14, 17, 13, 15};
Volume(2) = {2};
Surface Loop(3) = {18, 10, 7, 8, 9, 16};
Volume(3) = {3};

// ==== END [6um] ==== //
// Points and Lines
Point(23) = {0, 0, 6, 1.0};
// Ellipse 1
Point(24) = {0, -0.5, 6, 1.0};
Point(25) = {-0.5, -0.25, 6, 1.0};
Point(26) = {0.5, -0.25, 6, 1.0};
Point(27) = {0, -0.25, 6, 1.0};
// Ellipse 2
Point(28) = {0, 0.5, 6, 1.0};
Point(29) = {-0.5, 0.25, 6, 1.0};
Point(30) = {0.5, 0.25, 6, 1.0};
Point(31) = {0, 0.25, 6, 1.0};
// Major Axis Points
Point(32) = {0.1, -0.25, 6, 1.0};
Point(33) = {0.1, 0.25, 6, 1.0};
// Ellipse 1
Ellipse(30) = {25, 27, 32, 24};
Ellipse(31) = {24, 27, 32, 26};
Ellipse(32) = {26, 27, 32, 23};
Ellipse(33) = {23, 27, 32, 25};
// Ellipse 2
Ellipse(34) = {29, 31, 33, 23};
Ellipse(35) = {23, 31, 33, 30};
Ellipse(36) = {30, 31, 33, 28};
Ellipse(37) = {28, 31, 33, 29};
Line(38) = {14, 25};
Line(39) = {13, 24};
Line(40) = {15, 26};
Line(41) = {12, 23};
Line(42) = {18, 29};
Line(43) = {19, 30};
Line(44) = {17, 28};
// Loops and Surfaces
Curve Loop(55) = {30, 31, 32, 33};
Surface(21) = {55};
Curve Loop(57) = {34, 35, 36, 37};
Surface(22) = {57};

// ==== CONNECTION [4-6um] ==== //
Curve Loop(59) = {38, 30, -39, -13};
Surface(23) = {59};
Curve Loop(61) = {39, 31, -40, -14};
Surface(24) = {61};
Curve Loop(63) = {40, 32, -41, -15};
Surface(25) = {63};
Curve Loop(65) = {41, 33, -38, -16};
Surface(26) = {65};
Curve Loop(67) = {20, 42, -37, -44};
Surface(27) = {67};
Curve Loop(69) = {19, 44, -36, -43};
Surface(28) = {69};
Curve Loop(71) = {43, -35, -41, 18};
Surface(29) = {71};
Curve Loop(73) = {41, -34, -42, 17};
Surface(30) = {73};
Surface Loop(4) = {28, 27, 30, 22, 29, 17};
Volume(4) = {4};
Surface Loop(5) = {21, 23, 26, 24, 25, 18};
Volume(5) = {5};
