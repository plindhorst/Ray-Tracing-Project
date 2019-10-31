#include "BoundingBox.h"
#include "flyscene.hpp"

std::vector<BoundingBox*> BoundingBox::boxes = std::vector<BoundingBox*>();
Tucano::Mesh* BoundingBox::mesh = nullptr;

BoundingBox::BoundingBox(bool remember) {
	if (remember) boxes.push_back(this);
}

void BoundingBox::deconstruct() {
	for (BoundingBox* c : boxes) {
		delete c;
	}
	boxes.clear();
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
		if (!hasVertex(vertex)) { return false; }
	}
	return true;
}

bool BoundingBox::hasVertex(Eigen::Vector4f& vertex) {
	assert(vertex(3) == 1);
	return (vertex(0) >= low(0) && vertex(0) <= high(0)
		&& vertex(1) >= low(1) && vertex(1) <= high(1)
		&& vertex(2) >= low(2) && vertex(2) <= high(2));
}

void BoundingBox::fitMesh() {
	for (int i = 0; i < mesh->getNumberOfFaces(); i++) {
		faces.push_back(&mesh->getFace(i));
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
	faces.shrink_to_fit();
	outside.shrink_to_fit();
	return outside;
}

BoundingBox* BoundingBox::splitBox() {
	Eigen::Vector3f oldLow = low;
	Eigen::Vector3f oldHigh = high;
	int choice;

	if ((width >= height || failed[1]) && (width >= depth || failed[2]) && !failed[0]) {
		choice = 0;
	}
	else if ((height >= width || failed[0]) && (height >= depth || failed[2]) && !failed[1]) {
		choice = 1;
	}
	else if (!(failed[0] && failed[1] && failed[2])) {
		choice = 2;
	}
	else {
		return nullptr;
	}

	high(choice) = averageVertexCoord(choice);
	//high(choice) -= shape(choice) / 2;
	reshape();
	std::vector<Tucano::Face*> outsideFaces = this->outsideFaces();
	bool impossible = (faces.size() == 0 || outsideFaces.size() == 0);

	if (impossible) {
		if (faces.size() == 0) { faces = outsideFaces; }
		failed[choice] = true;
		low = oldLow;
		high = oldHigh;
		reshape();
		return this;
	}
	else {
		failed[0] = failed[1] = failed[2] = false;
		BoundingBox* newCube = new BoundingBox(true);
		newCube->faces = outsideFaces;
		newCube->fitFaces();
		fitFaces();
		return newCube;
	}
}

float BoundingBox::averageVertexCoord(int axis) {
	float average = 0;
	for (int i = 0; i < faces.size(); i++) {
		Tucano::Face* face = faces[i];
		average += mesh->getVertex(face->vertex_ids[0])(axis);
		average += mesh->getVertex(face->vertex_ids[1])(axis);
		average += mesh->getVertex(face->vertex_ids[2])(axis);
	}
	average /= faces.size() * 3;
	return average;
}

void BoundingBox::setRandomColor() {
	color = Eigen::Vector3f(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX);
}

int BoundingBox::getNumberOfFaces() {
	return faces.size();
}

void BoundingBox::generateShape() {
	box = new Tucano::Shapes::Box();
	Eigen::Affine3f modelMatrix = Eigen::Affine3f::Identity();
	Eigen::Affine3f meshMatrix = mesh->getShapeModelMatrix();
	Eigen::Vector3f tempLow = meshMatrix * low;
	Eigen::Vector3f tempHigh = meshMatrix * high;
	modelMatrix.translate(tempLow + 0.5 * (meshMatrix * high - meshMatrix * low));
	modelMatrix(0, 0) = tempHigh(0) - tempLow(0);
	modelMatrix(1, 1) = tempHigh(1) - tempLow(1);
	modelMatrix(2, 2) = tempHigh(2) - tempLow(2);
	box->setModelMatrix(modelMatrix);
	box->setColor(Eigen::Vector4f(color(0), color(1), color(2), 1));
}

void BoundingBox::render(Tucano::Flycamera& flyCamera, Tucano::Camera& scene_light) {
	if (faces.size() <= 0) { return; }
	box->render(flyCamera, scene_light);
}

BoundingBox::~BoundingBox() {
	delete box;
}