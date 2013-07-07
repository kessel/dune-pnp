
dna_refinement = 2.;
inner_pore_refinement = 2.;
outer_pore_refinement = 2.;
box_refinement = 6.;
pore_radius_left = 10;
pore_radius=10;

pore_length = 20;

dna_length = 20;

box_size_z = 100;
box_size_r = 55;

smoothing_radius = 1;

Point(0) = { -pore_length/2, box_size_r ,  0, outer_pore_refinement };
Point(1) = { -pore_length/2, pore_radius + smoothing_radius , 0, inner_pore_refinement };
Point(2) = { -pore_length/2+smoothing_radius, pore_radius , 0, inner_pore_refinement };
Point(3) = { -pore_length/2+smoothing_radius, pore_radius + smoothing_radius , 0, inner_pore_refinement };

//point(17) = { 0, 7, 0, 1 };

Point(4) = { +pore_length/2, box_size_r , 0, outer_pore_refinement };
Point(5) = { +pore_length/2, pore_radius_left + smoothing_radius , 0, inner_pore_refinement };
Point(6) = { +pore_length/2-smoothing_radius, pore_radius_left , 0, inner_pore_refinement };
Point(7) = { +pore_length/2-smoothing_radius, pore_radius_left + smoothing_radius , 0, inner_pore_refinement };

Line(111) = { 0 , 1 };
Line(1) = { 4 , 5 };
Line(2) = { 2, 6 };

Circle(3) = { 1, 3, 2 };
Circle(4) = { 5, 7, 6 };

Point(8) = { -dna_length/2, 0, 0, dna_refinement };
//Point(9) = { -dna_length/2 , dna_radius, 0, dna_refinement };
//Point(10) = { -dna_length/2 + dna_radius,  0, 0, dna_refinement };

Point(11) = { +dna_length/2, 0, 0, dna_refinement };
//Point(12) = { +dna_length/2 - dna_radius, dna_radius, 0, dna_refinement };
Point(13) = { +dna_length/2 ,  0, 0, dna_refinement };

//Circle(5) = { 8, 10, 9 };
//Circle(6) = { 11, 13, 12 };

Line(7) = { 8, 11 };

// Box
Point(14) = { - box_size_z/2, 0, 0, box_refinement };
Point(15) = { - box_size_z/2, box_size_r, 0, box_refinement };
Line(8) = { 8, 14 };
Line(9) = { 0, 15 };
Line(10) = { 14, 15 };

Point(16) = { + box_size_z/2, 0, 0, box_refinement };
Point(17) = { + box_size_z/2, box_size_r, 0, box_refinement };

Line(11) = { 11, 16 };
Line(12) = { 4, 17 };
Line(13) = { 16, 17 };




Line Loop(112) = {1, 4, -2, -3, -111, 9, -10, -8,  7,  11, 13, -12};
Plane Surface(113) = {112};
Physical Line(0) = {111, 3, 2, 4, 1};
Physical Line(1) = {11,7, 8};
Physical Line(2) = {10};
Physical Line(3) = {13};
Physical Line(4) = {9};
Physical Line(5) = {12};
Physical Surface(118) = {113};
