#include "mesh.h"
#include <fstream>

void split(string str, vector<string>& v, string spacer)
{
	int pos1, pos2;
	int len = spacer.length();     //记录分隔符的长度
	pos1 = 0;
	pos2 = str.find(spacer);
	while (pos2 != string::npos)
	{
		v.push_back(str.substr(pos1, pos2 - pos1));
		pos1 = pos2 + len;
		pos2 = str.find(spacer, pos1);    // 从str的pos1位置开始搜寻spacer
	}
	if (pos1 != str.length()) //分割最后一个部分
		v.push_back(str.substr(pos1));
}

void mesh2D::InitMesh(int type, double scale) {
	if (type == 0) {
		mesh_circle circle;
		for (int i = 0; i < circle.vertexNum; i++) {
			Vector2d vertex = scale * Vector2d(circle.vertex[i][0], circle.vertex[i][1]);
			Vector2d force = Vector2d(0, 0);
			Vector2d velocity = Vector2d(0, 0);
			double mass = 0;
			vertexes.push_back(vertex);
			forces.push_back(force);
			velocities.push_back(velocity);
			masses.push_back(mass);
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
			double mass = 0;
			vertexes.push_back(vertex);
			forces.push_back(force);
			velocities.push_back(velocity);
			masses.push_back(mass);
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
		Vector3d d_velocity = Vector3d(0, 0, 0);
		Vector3d d_pos = Vector3d(0, 0, 0);
		double mass = 0;
		Matrix3d Constraint; Constraint.setIdentity();
		vertexes.push_back(vertex);
		forces.push_back(force);
		velocities.push_back(velocity);
		d_velocities.push_back(d_velocity);
		Constraints.push_back(Constraint);
		masses.push_back(mass);
		d_positions.push_back(d_pos);
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


bool mesh3D::load_tetrahedraMesh(const std::string& filename, double scale) {

	ifstream ifs(filename);
	if (!ifs) {

		fprintf(stderr, "unable to read file %s\n", filename.c_str());
		ifs.close();
		exit(-1);
		return false;
	}

	double x, y, z;
	int index0, index1, index2, index3;
	string line = "";
	int nodeNumber = 0;
	int elementNumber = 0;
	while (getline(ifs, line)) {
		if (line.length() <= 1) continue;
		if (line == "$Nodes") {
			getline(ifs, line);
			nodeNumber = atoi(line.c_str());
			vertexNum = nodeNumber;
			for (int i = 0; i < nodeNumber; i++) {
				getline(ifs, line);
				vector<std::string> nodePos;
				std::string spacer = " ";
				split(line, nodePos, spacer);
				x = atof(nodePos[1].c_str());
				y = atof(nodePos[2].c_str());
				z = atof(nodePos[3].c_str());
				Vector3d d_velocity = Vector3d(0, 0, 0);
				Vector3d vertex = scale * Vector3d(x, y, z);
				Matrix3d Constraint; Constraint.setIdentity();
				Vector3d force = Vector3d(0, 0, 0);
				Vector3d velocity = Vector3d(0, 0, 0);
				Vector3d d_pos = Vector3d(0, 0, 0);
				double mass = 0;
				vertexes.push_back(vertex);
				forces.push_back(force);
				velocities.push_back(velocity);
				Constraints.push_back(Constraint);
				d_velocities.push_back(d_velocity);
				masses.push_back(mass);
				d_positions.push_back(d_pos);
			}
		}

		if (line == "$Elements") {
			getline(ifs, line);
			elementNumber = atoi(line.c_str());
			tetrahedraNum = elementNumber;
			for (int i = 0; i < elementNumber; i++) {
				getline(ifs, line);

				vector<std::string> elementIndexex;
				std::string spacer = " ";
				split(line, elementIndexex, spacer);
				index0 = atoi(elementIndexex[3].c_str()) - 1;
				index1 = atoi(elementIndexex[4].c_str()) - 1;
				index2 = atoi(elementIndexex[5].c_str()) - 1;
				index3 = atoi(elementIndexex[6].c_str()) - 1;

				vector<uint64_t> tetrahedra;
				tetrahedra.push_back(index0);
				tetrahedra.push_back(index1);
				tetrahedra.push_back(index2);
				tetrahedra.push_back(index3);
				tetrahedras.push_back(tetrahedra);
			}
			break;
		}
	}
	ifs.close();
	return true;
}

void generateCuboidMesh() {
	int length = 20, width = 20, height = 20;
	
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {

		}
	}
}