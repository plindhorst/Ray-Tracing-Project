#pragma once

#include <tucano/shapes/box.hpp>
#include <tucano/utils/flycamera.hpp>

class Box
{
public:
	static vector<Box*> boxes;
	Tucano::Shapes::Box mesh;

protected:
	Eigen::Vector3f low = Eigen::Vector3f(-0.5, -0.5, -0.5);
	Eigen::Vector3f high = Eigen::Vector3f(0.5, 0.5, 0.5);
	Eigen::Vector3f shape = high - low;
	float width = 1;
	float height = 1;
	float depth = 1;

	Eigen::Affine3f modelMatrix;
	vector<Tucano::Face*> faces;

public:
	Box();

	void fitMesh(Tucano::Mesh mesh);

	void render(Tucano::Flycamera& flycamera, Tucano::Camera& scene_light);

	void splitBox();

	void fit();

	void fit(Eigen::Vector3f low, Eigen::Vector3f high);

	// Renew Shape variables given low & high variables.
	void reshape();

protected:
	Box(const Box& otherBox);
};