#pragma once

#include <tucano/mesh.hpp>
#include <tucano/shapes/camerarep.hpp>
#include <tucano/shapes/cylinder.hpp>
#include <tucano/shapes/sphere.hpp>
#include <tucano/shapes/box.hpp>
#include <tucano/utils/flycamera.hpp>

class Cube
{
public:
	static std::vector<Cube*> cubes;
	Tucano::Shapes::Box box;

protected:
	Eigen::Vector3f low = Eigen::Vector3f(-0.5, -0.5, -0.5);
	Eigen::Vector3f high = Eigen::Vector3f(0.5, 0.5, 0.5);
	Eigen::Vector3f shape = high - low;
	float width = 1;
	float height = 1;
	float depth = 1;

	Tucano::Mesh* mesh;
	Eigen::Affine3f modelMatrix;
	vector<Tucano::Face*> faces;

public:
	Cube(bool remember);

	Cube(Eigen::Vector3f low, Eigen::Vector3f high, bool remember);

	void fitMesh(Tucano::Mesh& mesh);

	void render(const Tucano::Camera& flycamera, const Tucano::Camera& scene_light);

	void splitcube();

	void fit();

	void fit(Eigen::Vector3f low, Eigen::Vector3f high);

	// Renew Shape variables given low & high variables.
	void reshape();

	bool hasFace(Tucano::Face& face);

	bool hasVertex(Eigen::Vector4f& verted);

	// Deconstruct all created cubes.
	static void deconstruct();
};
