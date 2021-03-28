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
		mesh_rectangle rectangle;
		for (int i = 0; i < rectangle.vertexNum; i++) {
			Vector2d vertex = scale * Vector2d(rectangle.vertex[i][0], rectangle.vertex[i][1]);
			Vector2d force = Vector2d(0, 0);
			Vector2d velocity = Vector2d(0, 0);
			vertexes.push_back(vertex);
			forces.push_back(force);
			velocities.push_back(velocity);
		}

		for (int i = 0; i < rectangle.triangleNum; i++) {
			vector<uint64_t> triangle;
			triangle.push_back(rectangle.index[i][0]);
			triangle.push_back(rectangle.index[i][1]);
			triangle.push_back(rectangle.index[i][2]);
			//Vector3i triangle = Vector3i(circle.index[i][0], circle.index[i][1], circle.index[i][2]);
			triangles.push_back(triangle);
		}

		triangleNum = rectangle.triangleNum;
		vertexNum = rectangle.vertexNum;
	}
}

void mesh3D::InitMesh(int type, double scale) {
	mesh_cuboid cuboid;
	for (int i = 0; i < cuboid.vertexNum; i++) {
		Vector3d vertex = scale * Vector3d(cuboid.vertex[i][0], cuboid.vertex[i][1], cuboid.vertex[i][2]);
		Vector3d force = Vector3d(0, 0, 0);
		Vector3d velocity = Vector3d(0, 0, 0);
		vertexes.push_back(vertex);
		forces.push_back(force);
		velocities.push_back(velocity);
	}

	for (int i = 0; i < cuboid.tetrahedraNum; i++) {
		vector<uint64_t> tetrahedra;
		tetrahedra.push_back(cuboid.index[i][0]);
		tetrahedra.push_back(cuboid.index[i][1]);
		tetrahedra.push_back(cuboid.index[i][2]);
		tetrahedra.push_back(cuboid.index[i][3]);
		//Vector3i triangle = Vector3i(circle.index[i][0], circle.index[i][1], circle.index[i][2]);
		tetrahedras.push_back(tetrahedra);
	}

	tetrahedraNum = cuboid.tetrahedraNum;
	vertexNum = cuboid.vertexNum;
}