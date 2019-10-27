#include "BoundingBox.h"

std::vector<BoundingBox*> BoundingBox::boxes = std::vector<BoundingBox*>();
std::vector<Eigen::Vector3f> BoundingBox::triangleColors;
Tucano::Mesh* BoundingBox::mesh = nullptr;

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

void BoundingBox::fitMesh() {
	for (int i = 0; i < mesh->getNumberOfFaces(); i++) {
		faces.push_back(Face{ i, &mesh->getFace(i) });
	}
	fitFaces();
}

void BoundingBox::fitFaces() {
	// Predefine min and max coördinates to random vertex of random face.
	if (faces.size() <= 0) { return; }
	Eigen::Vector4f temp = mesh->getVertex(faces[0].face->vertex_ids[0]);
	float minx, miny, minz, maxx, maxy, maxz;
	minx = maxx = temp(0);
	miny = maxy = temp(1);
	minz = maxz = temp(2);

	// Find lowest and highest coördinates.
	for (int i = 0; i < faces.size(); i++) {
		Tucano::Face face = *faces[i].face;
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
	low = Eigen::Vector3f(minx, miny, minz);
	high = Eigen::Vector3f(maxx, maxy, maxz);
	reshape();
}

vector<Face> BoundingBox::outsideFaces() {
	vector<Face> outside = vector<Face>();
	vector<Face> inside = vector<Face>();
	for (auto i = faces.begin(); i != faces.end(); i++) {
		if (!hasFace(*(*i).face)) {
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
	if (width >= height && width >= depth) {
		//high(0) = averageVertexCoord(0);
		high(0) = high(0) - width / 2;
	}
	else if (height >= width && height >= depth) {
		//high(1) = averageVertexCoord(1);
		high(1) = high(1) - height / 2;
	}
	else {
		//high(2) = averageVertexCoord(2);
		high(2) = high(2) - depth / 2;
	}

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
	average /= mesh->getNumberOfVertices();
	return average;
}

void BoundingBox::setRandomColor(bool updateTriangleColors) {
	color = Eigen::Vector3f(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX);
	if (updateTriangleColors) {
		for (Face face : faces) {
			if (triangleColors.at(face.id)(0) == -1) {
				triangleColors.at(face.id) = color;
			}
			else {
				triangleColors.at(face.id) = (triangleColors.at(face.id) + color) / 2;
			}
		}
	}
}

int BoundingBox::getNumberOfFaces() {
	return faces.size();
}

void BoundingBox::generateShape() {
	box = Tucano::Shapes::Box();
	Eigen::Affine3f modelMatrix = Eigen::Affine3f::Identity();
	Eigen::Affine3f meshMatrix = mesh->getShapeModelMatrix();
	Eigen::Vector3f tempLow = meshMatrix * low;
	Eigen::Vector3f tempHigh = meshMatrix * high;
	modelMatrix.translate(tempLow + 0.5 * (meshMatrix * high - meshMatrix * low));
	modelMatrix(0, 0) = tempHigh(0) - tempLow(0);
	modelMatrix(1, 1) = tempHigh(1) - tempLow(1);
	modelMatrix(2, 2) = tempHigh(2) - tempLow(2);
	box.setModelMatrix(modelMatrix);
	box.setColor(Eigen::Vector4f(color(0), color(1), color(2), 1));
}

void BoundingBox::render(Tucano::Flycamera& flyCamera, Tucano::Camera& scene_light) {
	if (faces.size() <= 0) { return; }
	box.render(flyCamera, scene_light);
}