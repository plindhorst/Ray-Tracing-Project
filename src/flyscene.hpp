#ifndef __FLYSCENE__
#define __FLYSCENE__

// Must be included before glfw.
#include <GL/glew.h>

#include <GLFW/glfw3.h>

#include <tucano/effects/phongmaterialshader.hpp>
#include <tucano/mesh.hpp>
#include <tucano/shapes/camerarep.hpp>
#include <tucano/shapes/cylinder.hpp>
#include <tucano/shapes/sphere.hpp>
#include <tucano/shapes/box.hpp>
#include <tucano/utils/flycamera.hpp>
#include <tucano/utils/imageIO.hpp>
#include <tucano/utils/mtlIO.hpp>
#include <tucano/utils/objimporter.hpp>
#include <unordered_map>
#include "BoundingBox.h"

class Flyscene {

public:
	Flyscene(void) {}

	/**
	 * @brief Initializes the shader effect
	 * @param width Window width in pixels
	 * @param height Window height in pixels
	 */
	void initialize(int width, int height);

	/**
	 * Repaints screen buffer.
	 **/
	virtual void paintGL();

	/**
	 * Perform a single simulation step.
	 **/
	virtual void simulate(GLFWwindow* window);

	/**
	 * Returns the pointer to the flycamera instance
	 * @return pointer to flycamera
	 **/
	Tucano::Flycamera* getCamera(void) { return &flycamera; }

	/**
	 * @brief Add a new light source
	 */
	void addLight(void) { lights.push_back(flycamera.getCenter()); }

	/**
	 * @brief Create a debug ray at the current camera location and passing
	 * through pixel that mouse is over
	 * @param mouse_pos Mouse cursor position in pixels
	 */
	void createDebugRay(const Eigen::Vector2f& mouse_pos);

	/**
	 * @brief raytrace your scene from current camera position
	 */
	void raytraceScene(int width = 0, int height = 0);

	/**
	 * @brief trace a single ray from the camera passing through dest
	 * @param origin Ray origin
	 * @param dest Other point on the ray, usually screen coordinates
	 * @return a RGB color
	 */
	Eigen::Vector3f traceRay(Eigen::Vector3f& origin, Eigen::Vector3f& dir, int depth);

	//TO DO:: insert documentation
	void resetDebugRay();

	std::pair<Tucano::Face, std::pair<Eigen::Vector3f, float>> calculateMinimumFace(Eigen::Vector3f& origin, Eigen::Vector3f dir);

	// TO DO: insert documentation
	std::pair<Eigen::Vector3f, float> calculateDistance(Eigen::Vector3f& origin, Eigen::Vector3f& dir, Tucano::Face& face);

	// TO DO: insert documentation
	Eigen::Vector3f reflect(Eigen::Vector3f direction, Eigen::Vector3f normal);

	// TO DO: insert documentation
	Eigen::Vector3f refract(Eigen::Vector3f direction, Tucano::Face face);

	/**
	*Check if light is obstructed
	*/
	bool shadow(Eigen::Vector3f& dest, Eigen::Vector3f& light);

	void sphericalLight(Eigen::Vector3f& lightLoc, float radius, int nLightpoints);

	Eigen::Vector3f calcSingleColor(Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f& lightLoc, Eigen::Vector3f& pointP);

	Eigen::Vector3f calculateColor(Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f& pointP);

private:
	// A simple phong shader for rendering meshes
	Tucano::Effects::PhongMaterial phong;

	// A fly through camera
	Tucano::Flycamera flycamera;

	// the size of the image generated by ray tracing
	Eigen::Vector2i raytracing_image_size;

	// A camera representation for animating path (false means that we do not
	// render front face)
	Tucano::Shapes::CameraRep camerarep = Tucano::Shapes::CameraRep(false);

	// a frustum to represent the camera in the scene
	Tucano::Shapes::Sphere lightrep;

	// light sources for ray tracing
	vector<Eigen::Vector3f> lights;

	// Scene light represented as a camera
	Tucano::Camera scene_light;


	const static int max_depth = 2;

	/// A very thin cylinder to draw a debug ray (with the same value as max_depth)
	Tucano::Shapes::Cylinder ray [max_depth + 1];

	// Scene meshes
	Tucano::Mesh mesh;

	/// MTL materials
	vector<Tucano::Material::Mtl> materials;

	// SELFMADE

	void generateBoundingBoxes();
	void renderBoundingBoxes();

	bool intersectBox(Eigen::Vector3f& origin, Eigen::Vector3f& dir, BoundingBox& box);


public:
	static const bool RENDER_BOUNDINGBOXES = false;
	static const bool RENDER_BOUNDINGBOX_COLORED_TRIANGLES = false;
	static const int MIN_FACES = 300;

	const Eigen::Vector3f BACKGROUND_COLOR = Eigen::Vector3f(0.9, 0.9, 0.9);
	const Eigen::Vector3f FOREGROUND_COLOR = Eigen::Vector3f(0, 0, 1);

	const int nSphereLights = 30;

	// Default Material
	Eigen::Vector3f ka = Eigen::Vector3f(0.2, 0.2, 0.2);
	Eigen::Vector3f kd = Eigen::Vector3f(0.9, 0.9, 0);
	Eigen::Vector3f ks = Eigen::Vector3f(0, 0, 0);
	float shininess = 0;
	float refraction_index = 0;
	float transparency = 0;

	static std::unordered_map<Tucano::Face*, int> faceids;
	~Flyscene();
};

#endif // FLYSCENE
