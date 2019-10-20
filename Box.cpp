#include "Box.h"

Box::Box() {
	boxes.push_back(this);
	mesh = Tucano::Shapes::Box();
	modelMatrix = Eigen::Affine3f::Identity();
	mesh.setModelMatrix(modelMatrix);
}

void Box::fitMesh(Tucano::Mesh mesh) {
	float minx = FLT_MAX, miny = FLT_MAX, minz = FLT_MAX;
	float maxx = -FLT_MAX, maxy = -FLT_MAX, maxz = -FLT_MAX;

	for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
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
}

void Box::fit() {
	modelMatrix = Eigen::Affine3f::Identity();
	modelMatrix.translate(low + Eigen::Vector3f(0.5, 0.5, 0.5));
	modelMatrix(0, 0) = width;
	modelMatrix(1, 1) = height;
	modelMatrix(2, 2) = depth;
	mesh.setModelMatrix(modelMatrix);
}

void Box::fit(Eigen::Vector3f low, Eigen::Vector3f high) {
	this->low = low;
	this->high = high;
	reshape();
	fit();
}

void Box::reshape() {
	this->shape = high - low;
	width = shape(0);
	height = shape(1);
	depth = shape(2);
}

void Box::render(Tucano::Flycamera& flycamera, Tucano::Camera& scene_light) {
	mesh.render(flycamera, scene_light);
}

void Box::splitBox() {
	Eigen::Vector3f shift;

	if (width >= height && width >= depth) {
		width /= 2;
		shift = Eigen::Vector3f(width / 2, 0, 0);
	}
	else if (height >= width && height >= depth) {
		height /= 2;
		shift = Eigen::Vector3f(0, height / 2, 0);
	}
	else {
		depth /= 2;
		shift = Eigen::Vector3f(0, 0, depth / 2);
	}

	Box newBox = Box(*this);
	modelMatrix.translate(-shift);
	newBox.modelMatrix.translate(shift);
}