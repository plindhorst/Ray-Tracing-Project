#include "cube.h"

std::vector<Cube*> Cube::cubes = std::vector<Cube*>();

Cube::Cube(bool remember) {
	if (remember) cubes.push_back(this);
	box = Tucano::Shapes::Box();
	modelMatrix = Eigen::Affine3f::Identity();
	box.setModelMatrix(modelMatrix);
}

Cube::Cube(Eigen::Vector3f low, Eigen::Vector3f high, bool remember) {
	Cube(bool(remember));
	this->low = low;
	this->high = high;
	reshape();
}

void Cube::deconstruct() {
	for (Cube* c : cubes) {
		delete c;
	}
}

void Cube::render(const Tucano::Camera& flycamera, const Tucano::Camera& scene_light) {
	box.render(flycamera, scene_light);
}

void Cube::fit() {
	reshape();
	modelMatrix = Eigen::Affine3f::Identity();
	modelMatrix.translate(low + 0.5 * shape);
	modelMatrix(0, 0) = width;
	modelMatrix(1, 1) = height;
	modelMatrix(2, 2) = depth;
	box.setModelMatrix(modelMatrix);
}

void Cube::fit(Eigen::Vector3f low, Eigen::Vector3f high) {
	this->low = low;
	this->high = high;
	fit();
}

void Cube::reshape() {
	this->shape = high - low;
	width = shape(0);
	height = shape(1);
	depth = shape(2);
}

bool Cube::hasFace(Tucano::Face& face) {
	for (int id : face.vertex_ids) {
		Eigen::Vector4f vertex = mesh->getVertex(id);
		if (!hasVertex(vertex)) return false;
	}
	return true;
}

bool Cube::hasVertex(Eigen::Vector4f& vertex) {
	return (vertex(0) >= low(0) && vertex(0) <= high(0)
		&& vertex(1) >= low(1) && vertex(1) <= high(1)
		&& vertex(2) >= low(2) && vertex(2) <= high(2));
}

void Cube::fitMesh(Tucano::Mesh& mesh) {
	this->mesh = &mesh;
	for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
		faces.push_back(&mesh.getFace(i));
	}
	fitFaces();
}

void Cube::fitFaces() {
	// Predefine min and max coördinates to random vertex of random face.
	Eigen::Vector4f temp = mesh->getVertex(faces[0]->vertex_ids[0]);
	float minx, miny, minz, maxx, maxy, maxz;
	minx = maxx = temp(0);
	miny = maxy = temp(1);
	minz = maxz = temp(2);

	// Find lowest and highest coördinates.
	for (int i = 0; i < mesh->getNumberOfFaces(); i++) {
		Tucano::Face face = *faces[i];
		for (int j = 0; j < 3; j++) {
			Eigen::Vector4f v = mesh->getVertex(face.vertex_ids[j]);
			float x = v(0);
			float y = v(1);
			float z = v(2);
			if (x < minx) {
				minx = x;
			}
			else if (x > maxx) {
				maxx = x;
			}
			if (y < miny) {
				miny = y;
			}
			else if (y > maxy) {
				maxy = y;
			}
			if (z < minz) {
				minz = z;
			}
			else if (z > maxz) {
				maxz = z;
			}
		}
	}

	// Fit to coördinates.
	low = mesh->getShapeModelMatrix() * Eigen::Vector3f(minx, miny, minz);
	high = mesh->getShapeModelMatrix() * Eigen::Vector3f(maxx, maxy, maxz);
	fit();
}

void Cube::splitcube() {
	//cube* newcube;
	Eigen::Vector3f newHigh;

	if (width >= height && width >= depth) {
		width /= 2;
		newHigh = high - Eigen::Vector3f(-width, 0, 0);
	}
	else if (height >= width && height >= depth) {
		height /= 2;
		newHigh = high - Eigen::Vector3f(0, -height, 0);
	}
	else {
		depth /= 2;
		newHigh = high - Eigen::Vector3f(0, 0, -depth);
	}

	// Reshape first cube
	high = newHigh;
	reshape();
	fit();


}