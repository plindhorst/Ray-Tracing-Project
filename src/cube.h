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

protected:
	// Raw Coördinates (Object Space)
	Eigen::Vector3f rawLow;
	Eigen::Vector3f rawHigh;
	Eigen::Vector3f rawShape;
	float rawWidth;
	float rawHeight;
	float rawDepth;
	// Processed Coördinates (World Space)
	Eigen::Vector3f low = Eigen::Vector3f(-0.5, -0.5, -0.5);
	Eigen::Vector3f high = Eigen::Vector3f(0.5, 0.5, 0.5);
	Eigen::Vector3f shape = high - low;
	float width = 1;
	float height = 1;
	float depth = 1;

private:
	// The mesh this cube belongs to.
	Tucano::Mesh* mesh;
	// The box shape this cube represents.
	Tucano::Shapes::Box box;
	Eigen::Affine3f modelMatrix;
	// All faces that should be inside this cube.
	vector<Tucano::Face*> faces;

public:
	//	If remember, stores pointer to cube, will be deconstructed with class deconstructor.
	Cube(bool remember);

	// 	Delete all cubes that have been remembered.
	static void deconstruct();

	//	Render this cube using Tucano's box shape.
	void render(const Tucano::Camera& flycamera, const Tucano::Camera& scene_light);

	//	Sets ModelMatrix so the Cube fits its raw low & high coördinates.
	void fit();

	//	Renew shape variables given low & high variables.
	void reshape();

	//	Check if face is entirely inside this cube.
	bool hasFace(Tucano::Face& face);

	//	Check if this vertex is inside this cube (sides included).
	bool hasVertex(Eigen::Vector4f& verted);

	//	Fits Cube to contain entire mesh.
	void fitMesh(Tucano::Mesh& mesh);

	// 	Fits Cube to contain all it's faces.
	void fitFaces();

	//	Returns the faces that are outside the cube, and removes them from the list of faces that should be inside.
	vector<Tucano::Face*> outsideFaces();

	//	Split this box along longest axis, and create a new cube to cover all lost faces.
	void splitcube();

};
