	SetFactory("OpenCASCADE");
	Mesh.CharacteristicLengthMax = 10.0;
	Mesh.CharacteristicLengthMin = 10.0;
	R2 = 200;
	R1 = 100;

	Sphere(1) = {0,0,0, R2};
	Sphere(2) = {0,0,0, R1};
	BooleanDifference{ Volume{1}; Delete;}{ Volume{2}; Delete;}

	Physical Volume("hollow_volume", 10) = {1};
	Physical Surface("hollow_surface", 11) = {2, 3};
