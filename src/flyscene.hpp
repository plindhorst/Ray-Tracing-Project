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
	void addLight();

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

	void traceRayThread(int h, int w, int start, int stop, vector<vector<Eigen::Vector3f>>& pixel_data);

	//TO DO:: insert documentation
	void resetDebugRay();

	std::tuple<Tucano::Face, Eigen::Vector3f, float> calculateMinimumFace(Eigen::Vector3f& origin, Eigen::Vector3f dir);

	// TO DO: insert documentation
	std::pair<Eigen::Vector3f, float> calculateDistance(Eigen::Vector3f& origin, Eigen::Vector3f& dir, Tucano::Face& face);

	Eigen::Vector3f calculateReflectColor(Tucano::Face minimum_face, Eigen::Vector3f interPoint, Eigen::Vector3f& origin, Eigen::Vector3f& dir, int depth);

	// TO DO: insert documentation
	Eigen::Vector3f reflect(Eigen::Vector3f direction, Eigen::Vector3f normal);

	// TO DO: insert documentation
	Eigen::Vector3f refract(Eigen::Vector3f direction, Tucano::Face face);


	/**
	*Check if light is obstructed
	*/
	bool shadow(Eigen::Vector3f& dest, Eigen::Vector3f& light);

	void sphericalLight(std::pair<Eigen::Vector3f, Eigen::Vector3f> lightLoc, float radius, int nLightpoints);

	Eigen::Vector3f calcSingleColor(Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f lightDirection, Eigen::Vector3f light_intensity, Eigen::Vector3f& pointP);


	Eigen::Vector3f calculateColor(Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f& pointP);

private:
	int PIXEL_COUNT;

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

	Tucano::Shapes::Arrow dirLightrep;

	public: // so we can remove all lights
	// light sources for ray tracing
	vector<std::pair<Eigen::Vector3f, Eigen::Vector3f>> lights;

	// the directional lights
	vector<std::tuple<Eigen::Vector3f, Eigen::Vector3f, Eigen::Quaternion<float>>> dirLights;

	private:
	// Scene light represented as a camera
	Tucano::Camera scene_light;

	// Depth of recursive raytracing
	const static int max_depth = 2;

	//Array of cylinder for debugray
	Tucano::Shapes::Cylinder ray [max_depth];

	// Scene meshes
	Tucano::Mesh mesh;

	/// MTL materials
	vector<Tucano::Material::Mtl> materials;

	// SELFMADE

	void generateBoundingBoxes();
	void renderBoundingBoxes();

	bool intersectBox(Eigen::Vector3f& origin, Eigen::Vector3f& dir, BoundingBox& box);

	Eigen::Vector3f interpolateNormal(Tucano::Face& face, Eigen::Vector3f PointP);

public:
	const string OBJECT_NAME = "reflectivityy.obj";

	static const bool RENDER_BOUNDINGBOXES = false;
	static const bool RENDER_BOUNDINGBOX_COLORED_TRIANGLES = false;
	const int MIN_FACES = 300;
	const int MAX_BOXES = INT_MAX;

	static std::unordered_map<Tucano::Face*, int> faceids;

	static const int THREADS = 20;

	const Eigen::Vector3f BACKGROUND_COLOR = Eigen::Vector3f(0.9, 0.9, 0.9);
	const Eigen::Vector3f FOREGROUND_COLOR = Eigen::Vector3f(0, 0, 1);

	// Default Material
	Eigen::Vector3f ka = Eigen::Vector3f(0.2, 0.2, 0.2);
	Eigen::Vector3f kd = Eigen::Vector3f(0.9, 0.9, 0);
	Eigen::Vector3f ks = Eigen::Vector3f(0, 0, 0);
	float shininess = 0.f;
	float refraction_index = 0.f;
	float transparency = 0.f;
	~Flyscene();
};

#endif // FLYSCENE
