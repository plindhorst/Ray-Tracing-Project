#include "cube.h"

std::vector<Cube*> Cube::cubes = std::vector<Cube*>();

Cube::Cube(bool remember) {
	if (remember) cubes.push_back(this);
	box = Tucano::Shapes::Box();
	modelMatrix = Eigen::Affine3f::Identity();
	box.setModelMatrix(modelMatrix);
}

Cube::Cube(Eigen::Vector3f low, Eigen::Vector3f high, bool remember) {
	Cube(bool (remember));
	this->low = low;
	this->high = high;
	reshape();
}

void Cube::fitMesh(Tucano::Mesh& mesh) {
	Eigen::Vector4f temp = mesh.getVertex(0);
	float minx, miny, minz, maxx, maxy, maxz;
	minx = maxx = temp(0);
	miny = maxy = temp(1);
	minz = maxz = temp(2);

	for (int i = 1; i < mesh.getNumberOfVertices(); i++) {
		Eigen::Vector4f v = mesh.getVertex(i);
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

	low = mesh.getShapeModelMatrix() * Eigen::Vector3f(minx, miny, minz);
	high = mesh.getShapeModelMatrix() * Eigen::Vector3f(maxx, maxy, maxz);
	reshape();
	fit();

	for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
		faces.push_back(&mesh.getFace(i));
	}
	this->mesh = &mesh;
}

void Cube::fit() {
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
	reshape();
	fit();
}

void Cube::reshape() {
	this->shape = high - low;
	width = shape(0);
	height = shape(1);
	depth = shape(2);
}

void Cube::render(const Tucano::Camera& flycamera, const Tucano::Camera& scene_light) {
	box.render(flycamera, scene_light);
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

void Cube::deconstruct() {
	for (Cube* c : cubes) {
		delete c;
	}
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