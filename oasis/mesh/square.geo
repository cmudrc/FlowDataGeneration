Point(1) = {0, 0, 0, 0.0005};
Point(2) = {0.1, 0, 0, 0.0005};
Point(3) = {0.1, 0.07, 0, 0.0005};
Point(4) = {0., 0.07, 0, 0.0005};
Point(5) = {0.025, 0.03, 0, 0.00005};
Point(6) = {0.035, 0.03, 0, 0.00005};
Point(7) = {0.035, 0.04, 0, 0.00005};
Point(8) = {0.025, 0.04, 0, 0.00005};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 1};
Line(4) = {3, 4};
Line(5) = {8, 5};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {5, 6};
Line Loop(9) = {4, 3, 1, 2};
Line Loop(10) = {7, 5, 8, 6};
Plane Surface(11) = {9, 10};

