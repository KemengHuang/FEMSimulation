#pragma once
#ifndef FEM_MESH_H
#define FEM_MESH_H
#include "Eigen/Eigen"
#include <vector>
using namespace std;
using namespace Eigen;
struct mesh_circle {
	double vertex[9][2] = { {0,0},
		{1,0},
	{0.7071067811865476, 0.7071067811865476},
	{0, 1},
	{-0.7071067811865476, 0.7071067811865476},
	{-1, 0},
	{-0.7071067811865476, -0.7071067811865476},
	{0, -1},
	{0.7071067811865476, -0.7071067811865476} };

	//double force[9][2] = {0};
	//double velocity[9][2] = {0};

	int index[8][3] = { {1, 2, 0},
	{2, 3, 0},
	{3, 4, 0},
	{4, 5, 0},
	{5, 6, 0},
	{6, 7, 0},
	{7, 8, 0},
	{8, 1, 0} };
	//Matrix2d DM_triangle_inverse[8];
	int vertexNum = 9;
	int triangleNum = 8;
};

struct mesh_rectangle {
	double vertex[15][2] = { {0,1},
		{1,1},
		{1, 2},
		{0, 2},
		{-1, 2},
		{-1, 1},
		{-1, 0},
		{0, 0},
		{1, 0},
		{-1,-1},
		{0,-1},
		{1,-1},
		{-1,-2},
		{0, -2},
		{1,-2} };

	//double force[15][2] = { 0 };
	//double velocity[15][2] = { 0 };

	int index[16][3] = { {1, 2, 0},
		{2, 3, 0},
		{3, 4, 0},
		{4, 5, 0},
		{5, 6, 0},
		{6, 7, 0},
		{7, 8, 0},
		{8, 1, 0},
		{6, 7, 10},
		{7, 8, 10},
		{8, 11, 10},
		{11, 14, 10},
		{13, 14, 10},
		{12, 13, 10},
		{9, 12, 10},
		{6, 9, 10} };
	//Matrix2d DM_triangle_inverse[16];
	int vertexNum = 15;
	int triangleNum = 16;
};

struct mesh_cuboid {
	double vertex[8][3] = {
		{-1,1,1},
		{-1,1,-1},
		{1,1,-1},
		{1,1,1},
		{-1,-1,1},
		{-1,-1,-1},
		{1,-1,-1},
		{1,-1,1} };

	//double force[8][3] = { 0 };
	//double velocity[8][3] = { 0 };

	int index[5][4] = { {0, 4, 7, 5},
		{2, 3, 7, 0},
		{2, 6, 7, 5},
		{0, 1, 2, 5},
		{0, 2, 7, 5} };
	//Matrix3d DM_tetrahedra_inverse[5];
	int vertexNum = 8;
	int tetrahedraNum = 5;
};


class mesh2D {
public:
	vector<double> areas;
	vector<double> masses;
	vector<Vector2d> vertexes;
	vector<vector<uint64_t>> triangles;
	vector<Vector2d> forces;
	vector<Vector2d> velocities;
	vector<Matrix2d> DM_triangle_inverse;
	int vertexNum;
	int triangleNum;
	void InitMesh(int type, double scale);
};

class mesh3D {
public:
	double maxVolum;
	vector<Matrix3d> Constraints;
	vector<double> volum;
	vector<double> masses;
	vector<Vector3d> vertexes;
	vector<Vector3d> d_positions;
	vector<vector<uint64_t>> tetrahedras;
	vector<Vector3d> forces;
	vector<Vector3d> velocities;
	vector<Vector3d> d_velocities;
	vector<Matrix3d> DM_triangle_inverse;
	int vertexNum;
	int tetrahedraNum;
	void InitMesh(int type, double scale);
	bool load_tetrahedraMesh(const std::string& filename, double scale);
};

#endif // !FEM_MESH.H

