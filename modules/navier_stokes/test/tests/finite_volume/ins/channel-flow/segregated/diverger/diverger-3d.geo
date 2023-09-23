Point(1) = {0, 0, 0};
Point(2) = {0.5, 0, 0};
Point(3) = {0.7, -0.2, 0};
Point(4) = {1.3, -0.2, 0};
Point(5) = {1.3, 0.2, 0};
Point(6) = {0.7, 0.2, 0};
Point(7) = {0.5, 0.2, 0};
Point(8) = {0.0, 0.2, 0};

Line(1) = {1, 2}; Transfinite Curve{1} = 3;
Line(2) = {2, 3}; Transfinite Curve{2} = 3;
Line(3) = {3, 4}; Transfinite Curve{3} = 3;
Line(4) = {4, 5}; Transfinite Curve{4} = 3;
Line(5) = {5, 6}; Transfinite Curve{5} = 3;
Line(6) = {6, 7}; Transfinite Curve{6} = 3;
Line(7) = {7, 8}; Transfinite Curve{7} = 3;
Line(8) = {8, 1}; Transfinite Curve{8} = 3;
Line(9) = {2, 7}; Transfinite Curve{9} = 3;
Line(10) = {3, 6}; Transfinite Curve{10} = 3;

Curve Loop(1) = {1, 9, 7, 8};
Curve Loop(2) = {2, 10, 6, -9};
Curve Loop(3) = {3, 4, 5, -10};
Plane Surface(1) = {1}; Transfinite Surface{1};
Plane Surface(2) = {2}; Transfinite Surface{2};
Plane Surface(3) = {3}; Transfinite Surface{3};

new_surfs() = Extrude {0, 0, 0.1} { Surface{1, 2, 3}; Layers{1}; Recombine;};

Physical Volume("main") = {new_surfs[1], new_surfs[7], new_surfs[13]};

Physical Surface("front_and_back") = {1, 2, 3, new_surfs[0], new_surfs[6], new_surfs[12]};
Physical Surface("top") = {new_surfs[4], new_surfs[10], new_surfs[16]};
Physical Surface("bottom") = {new_surfs[2], new_surfs[8], new_surfs[14]};
Physical Surface("inlet") = {new_surfs[5]};
Physical Surface("outlet") = {new_surfs[15]};

Recombine Surface "*";
Recombine Volume "*";

