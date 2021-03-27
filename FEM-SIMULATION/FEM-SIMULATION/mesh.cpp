#include"mesh.h"
void mesh2D::InitMesh(int type, double scale) {
	if (type == 0) {
		mesh_circle circle;
		for (int i = 0; i < circle.vertexNum; i++) {
			Vector2d vertex = scale * Vector2d(circle.vertex[i][0], circle.vertex[i][1]);
			Vector2d force = Vector2d(0, 0);
			Vector2d velocity = Vector2d(0, 0);
			vertexes.push_back(vertex);
			forces.push_back(force);
			velocities.push_back(velocity);
		}

		for (int i = 0; i < circle.triangleNum; i++) {
			vector<uint64_t> triangle;
			triangle.push_back(circle.index[i][0]);
			triangle.push_back(circle.index[i][1]);
			triangle.push_back(circle.index[i][2]);
			//Vector3i triangle = Vector3i(circle.index[i][0], circle.index[i][1], circle.index[i][2]);
			triangles.push_back(triangle);
		}

		triangleNum = circle.triangleNum;
		vertexNum = circle.vertexNum;
	}
	else if (type == 1) {
		mesh_cuboid cuboid;
		for (int i = 0; i < cuboid.vertexNum; i++) {
			Vector2d vertex = scale * Vector2d(cuboid.vertex[i][0], cuboid.vertex[i][1]);
			Vector2d force = Vector2d(0, 0);
			Vector2d velocity = Vector2d(0, 0);
			vertexes.push_back(vertex);
			forces.push_back(force);
			velocities.push_back(velocity);
		}

		for (int i = 0; i < cuboid.triangleNum; i++) {
			vector<uint64_t> triangle;
			triangle.push_back(cuboid.index[i][0]);
			triangle.push_back(cuboid.index[i][1]);
			triangle.push_back(cuboid.index[i][2]);
			//Vector3i triangle = Vector3i(circle.index[i][0], circle.index[i][1], circle.index[i][2]);
			triangles.push_back(triangle);
		}

		triangleNum = cuboid.triangleNum;
		vertexNum = cuboid.vertexNum;
	}
}

void mesh3D::InitMesh(int type, double scale) {
	 
}