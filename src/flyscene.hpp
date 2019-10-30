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
	void addLight() {
		bool rightColor = false;
		float r, g, b;
		std::string l = "";
		while (!rightColor) {
			std::cout << " What colour do you want in rgb format? (Enter three numbers seperately between 1 and 0)" << std::endl;
			std::cin >> r >> g >> b;
			while (std::cin.fail()) {
				std::cout << "Error" << std::endl;
				std::cin.clear();
				std::cin.ignore(256, '\n');
				std::cin >> r >> g >> b;
			}
			if ((0 <= r <= 1) && (0 <= g <= 1) && (0 <= b <= 1)) {
				rightColor = true;
			}
			else {
				std::cout << " Wrong color code " << std::endl;
			}
		}
		while (!(l == "s" || l == "p" || l == "d") ) {
			std::cout << " Do you want a spherical, directional or a point light? (s, d or p)" << std::endl;
			std::cin >> l;
			while (std::cin.fail()) {
				std::cout << "Error" << std::endl;
				std::cin.clear();
				std::cin.ignore(256, '\n');
				std::cin >> l;
			}
		}
		if (l == "s" || l == "p") {
			std::pair<Eigen::Vector3f, Eigen::Vector3f> light = std::pair<Eigen::Vector3f, Eigen::Vector3f>(flycamera.getCenter(), Eigen::Vector3f(r, g, b));
			if (l == "s") {
				// create spherical lights out of point light
				string radius;
				std::cout << "What Radius?" << std::endl;
				std::cin >> radius;
				while (std::cin.fail()) {
					std::cout << "Error" << std::endl;
					std::cin.clear();
					std::cin.ignore(256, '\n');
					std::cin >> radius;
				}
				string lights;
				std::cout << "How many point light for approximations?" << std::endl;
				std::cin >> lights;
				while (std::cin.fail()) {
					std::cout << "Error" << std::endl;
					std::cin.clear();
					std::cin.ignore(256, '\n');
					std::cin >> lights;
				}
				int Nlights = std::stoi(lights);
				sphericalLight(light, std::stof(radius), Nlights);
				light.second /= (Nlights + 1);
			}

			lights.push_back(light);
		}

		if (l == "d") {
			std::cout << " What direction do you want the directional light to point at? " << std::endl;
			bool rightDir = false;
			float x, y, z;
			std::cin >> x >> y >> z;
			while (std::cin.fail()) {
				std::cout << "Error" << std::endl;
				std::cin.clear();
				std::cin.ignore(256, '\n');
				std::cin >> x >> y >> z;
			}
			std::pair<Eigen::Vector3f, Eigen::Vector3f> light = std::pair<Eigen::Vector3f, Eigen::Vector3f>(Eigen::Vector3f(x, y, z), Eigen::Vector3f(r, g, b));
		}
		std::cout << "Created a light with vector: (" << r << ", " << g << ", " << b << ")" << std::endl;
	}

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
	Eigen::Vector3f traceRay(Eigen::Vector3f& origin, Eigen::Vector3f& dir);

	void traceRayThread(int h, int w, int start, int stop, vector<vector<Eigen::Vector3f>>& pixel_data);

	// TO DO: insert documentation
	std::pair<Eigen::Vector3f, float> calculateDistance(Eigen::Vector3f& origin, Eigen::Vector3f& dir, Tucano::Face& face);

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

	public: // so we can remove all lights
	// light sources for ray tracing
	vector<std::pair<Eigen::Vector3f, Eigen::Vector3f>> lights;

	// the directional lights
	vector<std::pair<Eigen::Vector3f, Eigen::Vector3f>> dirLights;

	private:
	// Scene light represented as a camera
	Tucano::Camera scene_light;

	/// A very thin cylinder to draw a debug ray
	Tucano::Shapes::Cylinder ray = Tucano::Shapes::Cylinder(0.1, 1.0, 16, 64);

	// Scene meshes
	Tucano::Mesh mesh;

	/// MTL materials
	vector<Tucano::Material::Mtl> materials;

	// SELFMADE

	void generateBoundingBoxes();
	void renderBoundingBoxes();

	bool intersectBox(Eigen::Vector3f& origin, Eigen::Vector3f& dir, BoundingBox& box);

public:
	const string OBJECT_NAME = "dodgeColorTest.obj";

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
	float shininess = 0;
	float refraction_index = 0;
	float transparency = 0;
	~Flyscene();
};

#endif // FLYSCENE
