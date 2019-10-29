#pragma once

#include <tucano/mesh.hpp>
#include <tucano/shapes/camerarep.hpp>
#include <tucano/shapes/cylinder.hpp>
#include <tucano/shapes/sphere.hpp>
#include <tucano/shapes/box.hpp>
#include <tucano/utils/flycamera.hpp>

class BoundingBox
{
	public:
	static std::vector<BoundingBox*> boxes;
	static std::vector<Eigen::Vector3f> triangleColors;
	static Tucano::Mesh* mesh;

	// Raw Coördinates (Object Space)
	Eigen::Vector3f low;
	Eigen::Vector3f high;
	Eigen::Vector3f shape;
	float width;
	float height;
	float depth;
	Eigen::Vector3f color = Eigen::Vector3f(1,1,1);
	Tucano::Shapes::Box* box = nullptr;
	// All faces that should be inside this cube.
	vector<Tucano::Face*> faces;

	//	If remember, stores pointer to cube, will be deconstructed with class deconstructor.
	BoundingBox(bool remember);

	~BoundingBox();

	// 	Delete all boxes that have been remembered.
	static void deconstruct();

	//	Renew shape variables given low & high variables.
	void reshape();

	//	Check if face is entirely inside this cube.
	bool hasFace(Tucano::Face& face);

	//	Check if this vertex is inside this cube (sides included).
	bool hasVertex(Eigen::Vector4f& verted);

	//	Fits BoundingBox to contain entire mesh.
	void fitMesh();

	// 	Fits BoundingBox to contain all it's faces.
	void fitFaces();

	//	Returns the faces that are outside the cube, and removes them from the list of faces that should be inside.
	vector<Tucano::Face*> outsideFaces();

	//	Split this box along longest axis, and create a new cube to cover all lost faces.
	BoundingBox& splitBox();

	float averageVertexCoord(int axis);

	void setRandomColor();

	int getNumberOfFaces();

	void generateShape();

	void render(Tucano::Flycamera& flyCamera, Tucano::Camera& scene_light);
};
