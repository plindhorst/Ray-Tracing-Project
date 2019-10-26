#include "BoundingBox.h"

std::vector<BoundingBox*> BoundingBox::boxes = std::vector<BoundingBox*>();

BoundingBox::BoundingBox(bool remember) {
	if (remember) boxes.push_back(this);
}

void BoundingBox::deconstruct() {
	for (BoundingBox* c : boxes) {
		delete c;
	}
}

void BoundingBox::reshape() {
	// Raw
	shape = high - low;
	width = shape(0);
	height = shape(1);
	depth = shape(2);
}

bool BoundingBox::hasFace(Tucano::Face& face) {
	for (int id : face.vertex_ids) {
		Eigen::Vector4f vertex = mesh->getVertex(id);
		if (!hasVertex(vertex)) return false;
	}
	return true;
}

bool BoundingBox::hasVertex(Eigen::Vector4f& vertex) {
	return (vertex(0) >= low(0) && vertex(0) <= high(0)
		&& vertex(1) >= low(1) && vertex(1) <= high(1)
		&& vertex(2) >= low(2) && vertex(2) <= high(2));
}

void BoundingBox::fitMesh(Tucano::Mesh& mesh) {
	this->mesh = &mesh;
	for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
		faces.push_back(&mesh.getFace(i));
	}
	fitFaces();
}

void BoundingBox::fitFaces() {
	// Predefine min and max co�rdinates to random vertex of random face.
	if (faces.size() <= 0) { return; }
	Eigen::Vector4f temp = mesh->getVertex(faces[0]->vertex_ids[0]);
	float minx, miny, minz, maxx, maxy, maxz;
	minx = maxx = temp(0);
	miny = maxy = temp(1);
	minz = maxz = temp(2);

	// Find lowest and highest co�rdinates.
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

	// Fit to co�rdinates.
	low = Eigen::Vector3f(minx, miny, minz);
	high = Eigen::Vector3f(maxx, maxy, maxz);
	reshape();
}

vector<Tucano::Face*> BoundingBox::outsideFaces() {
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

BoundingBox& BoundingBox::splitBox() {
	Eigen::Vector3f newHigh;

	if (width >= height && width >= depth) {
		high(0) = averageVertexCoord(0);
	}
	else if (height >= width && height >= depth) {
		high(1) = averageVertexCoord(1);
	}
	else {
		high(2) = averageVertexCoord(2);
	}

	// Reshape first cube
	reshape();

	BoundingBox* newCube = new BoundingBox(true);
	newCube->mesh = mesh;
	newCube->faces = outsideFaces();
	newCube->fitFaces();
	fitFaces();
	return *newCube;
}

float BoundingBox::averageVertexCoord(int axis) {
	float average = 0;
	for (int i = 0; i < mesh->getNumberOfVertices(); i++) {
		average += mesh->getVertex(i)(axis);
	}
	average /= mesh->getNumberOfVertices;
	return average;
}

void BoundingBox::setRandomColor() {
	color = Eigen::Vector4f(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX);
}

int BoundingBox::getNumberOfFaces() {
	return faces.size();
}