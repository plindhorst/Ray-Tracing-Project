#include "cube.h"

std::vector<Cube*> Cube::cubes = std::vector<Cube*>();

Cube::Cube(bool remember) {
	if (remember) cubes.push_back(this);
	box = Tucano::Shapes::Box();
	modelMatrix = Eigen::Affine3f::Identity();
	box.setModelMatrix(modelMatrix);
}

void Cube::deconstruct() {
	for (Cube* c : cubes) {
		delete c;
	}
}

void Cube::render(const Tucano::Camera& flycamera, const Tucano::Camera& scene_light) {
	if (faces.size() <= 0) { return; }
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

void Cube::reshape() {
	// Raw
	rawShape = rawHigh - rawLow;
	rawWidth = rawShape(0);
	rawHeight = rawShape(1);
	rawDepth = rawShape(2);
	// Processed
	low = mesh->getShapeModelMatrix() * rawLow;
	high = mesh->getShapeModelMatrix() * rawHigh;
	shape = high - low;
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
	return (vertex(0) >= rawLow(0) && vertex(0) <= rawHigh(0)
		&& vertex(1) >= rawLow(1) && vertex(1) <= rawHigh(1)
		&& vertex(2) >= rawLow(2) && vertex(2) <= rawHigh(2));
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
	if (faces.size() <= 0) { return; }
	Eigen::Vector4f temp = mesh->getVertex(faces[0]->vertex_ids[0]);
	float minx, miny, minz, maxx, maxy, maxz;
	minx = maxx = temp(0);
	miny = maxy = temp(1);
	minz = maxz = temp(2);

	// Find lowest and highest coördinates.
	for (int i = 0; i < faces.size(); i++) {
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
	rawLow = Eigen::Vector3f(minx, miny, minz);
	rawHigh = Eigen::Vector3f(maxx, maxy, maxz);
	fit();
}

vector<Tucano::Face*> Cube::outsideFaces() {
	vector<Tucano::Face*> outside = vector<Tucano::Face*>();
	vector<Tucano::Face*> inside = vector<Tucano::Face*>();
	for (auto i = faces.begin(); i != faces.end(); i++) {
		if (!hasFace(**i)) {
			outside.push_back(*i);
		}
		else {
			inside.push_back(*i);
		}
	}
	faces = inside;
	return outside;
}

Cube& Cube::splitcube() {
	Eigen::Vector3f newRawHigh;

	if (width >= height && width >= depth) {
		rawWidth /= 2;
		newRawHigh = rawHigh - Eigen::Vector3f(rawWidth, 0, 0);
	}
	else if (height >= width && height >= depth) {
		rawHeight /= 2;
		newRawHigh = rawHigh - Eigen::Vector3f(0, rawHeight, 0);
	}
	else {
		rawDepth /= 2;
		newRawHigh = rawHigh - Eigen::Vector3f(0, 0, rawDepth);
	}

	// Reshape first cube
	rawHigh = newRawHigh;
	reshape();

	Cube* newCube = new Cube(true);
	newCube->mesh = mesh;
	newCube->faces = outsideFaces();
	newCube->fitFaces();
	newCube->box.setColor(Eigen::Vector4f(1, 0, 0, 0.5));
	fitFaces();
	return *newCube;
}

Eigen::Vector4f Cube::setRandomColor() {
	Eigen::Vector4f color = Eigen::Vector4f(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, 1);
	box.setColor(color);
	return color;
}

int Cube::getNumberOfFaces() {
	return faces.size();
}