#include <algorithm>
#include <thread>
#include <time.h>

#include "BoundingBox.h"
#include "flyscene.hpp"
#include <GLFW/glfw3.h>

std::unordered_map<Tucano::Face*, int> Flyscene::faceids;

void Flyscene::initialize(int width, int height) {
	// initiliaze the Phong Shading effect for the Opengl Previewer
	phong.initialize();

	// set the camera's projection matrix
	flycamera.setPerspectiveMatrix(60.0, width / (float)height, 0.1f, 100.0f);
	flycamera.setViewport(Eigen::Vector2f((float)width, (float)height));

	// load the OBJ file and materials
	Tucano::MeshImporter::loadObjFile(mesh, materials,
		"resources/models/" + OBJECT_NAME);


	// normalize the model (scale to unit cube and center at origin)
	mesh.normalizeModelMatrix();

	// pass all the materials to the Phong Shader
	for (int i = 0; i < materials.size(); ++i)
		phong.addMaterial(materials[i]);

	// set the color and size of the sphere to represent the light sources
	// same sphere is used for all sources
	lightrep.setColor(Eigen::Vector4f(1.0, 1.0, 0.0, 1.0));
	lightrep.setSize(0.15);

	dirLightrep.setColor(Eigen::Vector4f(1, 1, 0, 1));
	dirLightrep.setSize(0.15);

	// create a first ray-tracing light source at some random position
	lights.push_back(std::pair<Eigen::Vector3f, Eigen::Vector3f>(Eigen::Vector3f(-0.5, 2.0, 3.0), Eigen::Vector3f(0, 0, 0)));

	// scale the camera representation (frustum) for the ray debug
	camerarep.shapeMatrix()->scale(0.2);

	// the debug ray is a cylinder, set the radius and length of the cylinder
	ray.setSize(0.005, 10.0);

	// craete a first debug ray pointing at the center of the screen
	createDebugRay(Eigen::Vector2f(width / 2.0, height / 2.0));

	glEnable(GL_DEPTH_TEST);

	// for (int i = 0; i<mesh.getNumberOfFaces(); ++i){
	//   Tucano::Face face = mesh.getFace(i);    
	//   for (int j =0; j<face.vertex_ids.size(); ++j){
	//     std::cout<<"vid "<<j<<" "<<face.vertex_ids[j]<<std::endl;
	//     std::cout<<"vertex "<<mesh.getVertex(face.vertex_ids[j]).transpose()<<std::endl;
	//     std::cout<<"normal "<<mesh.getNormal(face.vertex_ids[j]).transpose()<<std::endl;
	//   }
	//   std::cout<<"mat id "<<face.material_id<<std::endl<<std::endl;
	//   std::cout<<"face   normal "<<face.normal.transpose() << std::endl << std::endl;
	// }

	generateBoundingBoxes();

	if (RENDER_BOUNDINGBOX_COLORED_TRIANGLES) {
		// Create face pointer to id map.
		faceids = std::unordered_map<Tucano::Face*, int>(mesh.getNumberOfFaces());
		for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
			faceids.insert(std::pair<Tucano::Face*, int>(&mesh.getFace(i), i));
		}
	}
}

void Flyscene::paintGL(void) {

	// update the camera view matrix with the last mouse interactions
	flycamera.updateViewMatrix();
	Eigen::Vector4f viewport = flycamera.getViewport();

	// clear the screen and set background color
	glClearColor(0.9, 0.9, 0.9, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// position the scene light at the last ray-tracing light source
	scene_light.resetViewMatrix();
	scene_light.viewMatrix()->translate(-lights.back().first);

	// render the scene using OpenGL and one light source
	phong.render(mesh, flycamera, scene_light);

	// render the ray and camera representation for ray debug
	ray.render(flycamera, scene_light);
	camerarep.render(flycamera, scene_light);

	// render ray tracing light sources as yellow spheres
	for (int i = 0; i < lights.size(); ++i) {
		lightrep.resetModelMatrix();
		lightrep.modelMatrix()->translate(lights[i].first);
		lightrep.render(flycamera, scene_light);
	}
	for (int i = 0; i < dirLights.size(); ++i) {

		dirLightrep.resetModelMatrix();
		Eigen::Affine3f m = dirLightrep.getModelMatrix();
		m(2, 2) = -1;
		dirLightrep.setModelMatrix(m);
		Eigen::Quaternion<float> rotation = get<2>(dirLights[i]);
		dirLightrep.modelMatrix()->rotate(rotation);
		dirLightrep.modelMatrix()->translate(Eigen::Vector3f(0, 1, 0));
		dirLightrep.render(flycamera, scene_light);
	}

	// render coordinate system at lower right corner
	flycamera.renderAtCorner();

	renderBoundingBoxes();
}

void Flyscene::simulate(GLFWwindow* window) {
	// Update the camera.
	// NOTE(mickvangelderen): GLFW 3.2 has a problem on ubuntu where some key
	// events are repeated: https://github.com/glfw/glfw/issues/747. Sucks.
	float dx = (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS ? 0.5 : 0.0) -
		(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS ? 0.5 : 0.0);
	float dy = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ? 0.5 : 0.0) -
		(glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ? 0.5 : 0.0);
	float dz = (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS ? 0.5 : 0.0) -
		(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS ? 0.5 : 0.0);
	flycamera.translate(dx, dy, dz);
}

void Flyscene::createDebugRay(const Eigen::Vector2f& mouse_pos) {
	ray.resetModelMatrix();
	// from pixel position to world coordinates
	Eigen::Vector3f screen_pos = flycamera.screenToWorld(mouse_pos);

	// direction from camera center to click position
	Eigen::Vector3f dir = (screen_pos - flycamera.getCenter()).normalized();

	// position and orient the cylinder representing the ray
	ray.setOriginOrientation(flycamera.getCenter(), dir);

	// place the camera representation (frustum) on current camera location, 
	camerarep.resetModelMatrix();
	camerarep.setModelMatrix(flycamera.getViewMatrix().inverse());
}

void Flyscene::addLight() {
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
	while (!(l == "s" || l == "p" || l == "d")) {
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
		Eigen::Vector3f direction = flycamera.screenToWorld(Eigen::Vector2f(flycamera.getViewport()(2) / 2, flycamera.getViewport()(3) / 2));
		std::tuple<Eigen::Vector3f, Eigen::Vector3f, Eigen::Quaternion<float>> light = std::tuple<Eigen::Vector3f, Eigen::Vector3f, Eigen::Quaternion<float>>(direction, Eigen::Vector3f(r, g, b), flycamera.getRotationMatrix());
		dirLights.push_back(light);
	}
	std::cout << "Created a light with color spectrum (r,g,b): (" << r << ", " << g << ", " << b << ")" << std::endl;
}

void Flyscene::raytraceScene(int width, int height) {
	std::cout << "ray tracing ..." << std::endl;
	time_t begin_time = time(NULL);

	// if no width or height passed, use dimensions of current viewport
	Eigen::Vector2i image_size(width, height);
	if (width == 0 || height == 0) {
		image_size = flycamera.getViewportSize();
	}

	// create 2d vector to hold pixel colors and resize to match image size
	vector<vector<Eigen::Vector3f>> pixel_data;
	pixel_data.resize(image_size[1]);
	for (int i = 0; i < image_size[1]; ++i)
		pixel_data[i].resize(image_size[0]);

	std::vector<std::thread> threads;

	// Check if image size can be divided by the number of threads
	if (image_size[1] % THREADS != 0) {
		std::cout << "Error: Incorrect number of threads" << std::endl;
		return;
	}

	PIXEL_COUNT = 0;
	int start = 0;
	int pixels = image_size[1] / THREADS;
	for (int i = 0; i < THREADS; i++) {
		threads.emplace_back(&Flyscene::traceRayThread, this, image_size[1], image_size[0], start, start + pixels, std::ref(pixel_data));
		start += pixels;
	}

	while (PIXEL_COUNT < 1000) {
		std::cout << "\r" << (PIXEL_COUNT + 1) << "/" << image_size[1];
	}

	// Wait for all threads to end
	for (std::thread& t : threads) {
		t.join();
	}

	// write the ray tracing result to a PPM image
	Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);

	time_t end_time = time(NULL);
	std::cout << " ray tracing done! " << std::endl;
	std::cout << std::endl << "Time it took to render(in seconds): " << end_time - begin_time << std::endl;
}

void Flyscene::traceRayThread(int h, int w, int start, int stop, vector<vector<Eigen::Vector3f>>& pixel_data) {
	// origin of the ray is always the camera center
	Eigen::Vector3f origin = flycamera.getCenter();
	Eigen::Vector3f direction;

	// for every pixel shoot a ray from the origin through the pixel coords
	for (int j = start; j < stop; ++j) {
		for (int i = 0; i < h; ++i) {
			// create a ray from the camera passing through the pixel (i,j)
			direction = (flycamera.screenToWorld(Eigen::Vector2f(i, j)) - origin).normalized();
			// launch raytracing for the given ray and write result to pixel data
			pixel_data[j][i] = traceRay(origin, direction);
		}
		PIXEL_COUNT++;
	}
}


Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f& origin, Eigen::Vector3f& dir) {
	// Parameters to keep track of current faces and the closest face
	std::pair<Eigen::Vector3f, float> minimum_distance_and_point = std::pair<Eigen::Vector3f, float>(Eigen::Vector3f::Zero(), INFINITY);
	Tucano::Face minimum_face;
	std::pair<Eigen::Vector3f, float> current_distance_and_point = std::pair<Eigen::Vector3f, float>(Eigen::Vector3f::Zero(), INFINITY);
	Tucano::Face current_face;

	// Loop through all Bounding boxes.
	for (BoundingBox* box : BoundingBox::boxes) {
		if (intersectBox(origin, dir, *box)) {
			for (int i = 0; i < box->faces.size(); i++) {
				current_face = *box->faces[i];
				current_distance_and_point = calculateDistance(origin, dir, current_face);
				if (0 <= current_distance_and_point.second && current_distance_and_point.second < minimum_distance_and_point.second) {
					minimum_distance_and_point = current_distance_and_point;
					minimum_face = current_face;
				}
			}
		}
	}
	// Test if the ray intersected with a face, if so: calculate the color
	if (minimum_distance_and_point.second == INFINITY) {
		return BACKGROUND_COLOR;
	}
	if (RENDER_BOUNDINGBOX_COLORED_TRIANGLES) {
		return BoundingBox::triangleColors.at(faceids[&minimum_face]);
	}
	return calculateColor(minimum_face, origin, minimum_distance_and_point.first);
}

void Flyscene::generateBoundingBoxes() {
	BoundingBox::mesh = &mesh;
	BoundingBox::triangleColors = std::vector<Eigen::Vector3f>(mesh.getNumberOfFaces(), Eigen::Vector3f(-1, -1, -1));
	BoundingBox* box = new BoundingBox(true);
	box->fitMesh();

	bool notDone = true;
	while (notDone && BoundingBox::boxes.size() < MAX_BOXES) {
		notDone = false;
		vector<BoundingBox*> current = BoundingBox::boxes;
		std::cout << "Now have " << current.size() << " boxes." << std::endl;
		for (BoundingBox* box : current) {
			if (box->getNumberOfFaces() > MIN_FACES && (!box->failed[0] || !box->failed[1] || !box->failed[2])) {
				BoundingBox* newBox = box->splitBox();
				while (box == newBox) {
					newBox = box->splitBox();
				}
				notDone = true;
			}
		}
	}

	std::cout << "Boxes: " << BoundingBox::boxes.size() << std::endl;

	if (RENDER_BOUNDINGBOXES || RENDER_BOUNDINGBOX_COLORED_TRIANGLES) {
		for (BoundingBox* box : BoundingBox::boxes) {
			box->setRandomColor();
			box->generateShape();
		}
	}
}

void Flyscene::renderBoundingBoxes() {
	if (!RENDER_BOUNDINGBOXES) {
		return;
	}
	for (BoundingBox* boundingBox : BoundingBox::boxes) {
		boundingBox->render(flycamera, scene_light);
	}
}

Flyscene::~Flyscene() {
	BoundingBox::deconstruct();
}


std::pair<Eigen::Vector3f, float> Flyscene::calculateDistance(Eigen::Vector3f& origin, Eigen::Vector3f& dir, Tucano::Face& face) {
	//get a vertex and the normal of the face
	Eigen::Affine3f modelMatrix = mesh.getShapeModelMatrix();

	Eigen::Vector4f vec0 = mesh.getVertex(face.vertex_ids[0]);
	Eigen::Vector3f vert0 = modelMatrix * (vec0.head<3>() / vec0.w()); // Note: w should always be 1, which may be assumed.
	Eigen::Vector3f facenormal = face.normal.normalized();

	//Return false if triangle and direction of ray are the same
	if (facenormal.dot(dir) == 0) {
		std::pair<Eigen::Vector3f, float> pair;
		pair.second = (float)-1;
		return pair;
	}
	//Calculate the distance between the plane and the origin (not the camera)
	float distancePlane = facenormal.dot(vert0);
	//Ray = origin + t*distance

	float orthProjectionDest = distancePlane - origin.dot(facenormal);
	float t = orthProjectionDest / (dir.dot(facenormal));
	Eigen::Vector3f PointP = origin + t * dir;

	//If any area is smaller or equal to zero, point is on the wrong side of the edge. 
	if (interpolateNormal(face, PointP).norm() == 0) {
		std::pair<Eigen::Vector3f, float> pair;
		pair.second = (float)-1;
		return pair;
	}
	else {
		std::pair<Eigen::Vector3f, float> pair;
		pair.first = PointP;
		pair.second = t;
		return pair;
	}
}

bool Flyscene::intersectBox(Eigen::Vector3f& origin, Eigen::Vector3f& dir, BoundingBox& box) {
	Eigen::Affine3f ME = mesh.getShapeModelMatrix();
	Eigen::Affine3f M = mesh.getShapeModelMatrix().inverse();
	Eigen::Matrix3f MS = Eigen::Matrix3f(M.matrix().block<3, 3>(0, 0));

	Eigen::Vector3f origin2 = M * origin;
	Eigen::Vector3f dir2 = (MS * dir).normalized();
	Eigen::Vector3f tmin = box.low - origin2;
	Eigen::Vector3f tmax = box.high - origin2;
	tmin(0) /= dir2(0);
	tmax(0) /= dir2(0);
	tmin(1) /= dir2(1);
	tmax(1) /= dir2(1);
	tmin(2) /= dir2(2);
	tmax(2) /= dir2(2);

	Eigen::Vector3f tinv = Eigen::Vector3f(min(tmin(0), tmax(0)), min(tmin(1), tmax(1)), min(tmin(2), tmax(2)));
	Eigen::Vector3f toutv = Eigen::Vector3f(max(tmin(0), tmax(0)), max(tmin(1), tmax(1)), max(tmin(2), tmax(2)));

	float tin = max(tinv(0), max(tinv(1), tinv(2)));
	float tout = min(toutv(0), min(toutv(1), toutv(2)));

	return !(tin > tout || tout < 0);
}

//Call function in calculate color, if it returns true => pixel should be black/ambient. If it returns false => the pixel should have a color.
bool Flyscene::shadow(Eigen::Vector3f& pointP, Eigen::Vector3f& lightDirection) {
	Tucano::Face current_face;
	Eigen::Vector3f inter = pointP + 0.001 * (lightDirection);

	// Loop through all Bounding boxes.
	for (BoundingBox* box : BoundingBox::boxes) {
		if (intersectBox(pointP, lightDirection, *box)) {
			for (int i = 0; i < box->faces.size(); i++) {
				current_face = *box->faces[i];
				if (calculateDistance(inter, lightDirection, current_face).second >= 0) {
					return true;
				}
			}
		}
	}
	return false;
}

//Call function in rayTraceScene before pixel loop
void Flyscene::sphericalLight(std::pair<Eigen::Vector3f, Eigen::Vector3f> light, float radius, int nLightpoints) {
	Eigen::Vector3f lightLoc = light.first;
	for (int i = 0; i < nLightpoints; i++) {
		float x = (-radius + (rand() / (RAND_MAX / (radius * 2))));
		float y = (-radius + (rand() / (RAND_MAX / (radius * 2))));
		float z = (-radius + (rand() / (RAND_MAX / (radius * 2))));
		lights.push_back(std::pair<Eigen::Vector3f, Eigen::Vector3f>(Eigen::Vector3f(lightLoc.x() + x, lightLoc.y() + y, lightLoc.z() + z), light.second / (nLightpoints + 1)));
	}
}


//Call function in calculateColor
Eigen::Vector3f Flyscene::calcSingleColor(Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f lightDirection, Eigen::Vector3f light_intensity, Eigen::Vector3f& pointP) {
	if (shadow(pointP, lightDirection)) {
		return Eigen::Vector3f(0.0, 0.0, 0.0);
	}
	if (minimum_face.material_id != -1) {
		Tucano::Material::Mtl material = materials[minimum_face.material_id];
		ka = material.getAmbient();
		kd = material.getDiffuse();
		ks = material.getSpecular();
		shininess = material.getShininess();
		refraction_index = material.getOpticalDensity();
		transparency = material.getDissolveFactor();
	}

	Eigen::Vector3f normal = interpolateNormal(minimum_face, PointP);
	Eigen::Vector3f lightReflection = (lightDirection - 2 * (normal.dot(lightDirection)) * normal);
	Eigen::Vector3f eyeDirection = (origin - pointP).normalized();

	Eigen::Vector3f ambient = Eigen::Vector3f(light_intensity.x() * ka.x(), light_intensity.y() * ka.y(), light_intensity.z() * ka.z());
	Eigen::Vector3f diffuse = Eigen::Vector3f(light_intensity.x() * kd.x(), light_intensity.y() * kd.y(), light_intensity.z() * kd.z()) * std::max(lightDirection.dot(normal), 0.f);
	Eigen::Vector3f specular = Eigen::Vector3f(light_intensity.x() * ks.x(), light_intensity.y() * ks.y(), light_intensity.z() * ks.z()) * std::max(std::pow(lightReflection.dot(eyeDirection), shininess), 0.f);

	Eigen::Vector3f color = ambient + diffuse + specular;
	return color;
}


/**
	Returns an interpolated normal of the point P on the face. Returns nullvector if point does not lie on the face.
*/
Eigen::Vector3f Flyscene::interpolateNormal(Tucano::Face& face, Eigen::Vector3f PointP) {
	Eigen::Vector3f vert0 = mesh.getVertex(face.vertex_ids[0]);
	Eigen::Vector3f vert1 = mesh.getVertex(face.vertex_ids[1]);
	Eigen::Vector3f vert2 = mesh.getVertex(face.vertex_ids[2]);
	Eigen::Vector3f facenormal = face.normal.normalized();

	Eigen::Vector3f edge0 = vert1 - vert0;
	Eigen::Vector3f edge1 = vert2 - vert1;
	Eigen::Vector3f edge2 = vert0 - vert2;

	Eigen::Vector3f Inner0 = PointP - vert0;
	Eigen::Vector3f Inner1 = PointP - vert1;
	Eigen::Vector3f Inner2 = PointP - vert2;

	float area0 = edge0.cross(Inner0).norm() / 2;
	float area1 = edge1.cross(Inner1).norm() / 2;
	float area2 = edge2.cross(Inner2).norm() / 2;

	if (area0 < 0 || area1 < 0 || area2 < 0) {
		return Eigen::Vector3f(0, 0, 0);
	}

	float area = edge0.cross(-edge2).norm() / 2;
	return (mesh.getNormal(face.vertex_ids[0]).normalized() * area1 / area + mesh.getNormal(face.vertex_ids[1]).normalized() * area2 / area + mesh.getNormal(face.vertex_ids[2]).normalized() * area0 / area).normalized();
}

//Call function in traceRay
Eigen::Vector3f Flyscene::calculateColor(Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f& pointP) {
	Eigen::Vector3f sumColor = Eigen::Vector3f(0.0, 0.0, 0.0);
	Eigen::Vector3f sumColorDir = Eigen::Vector3f(0.0, 0.0, 0.0);
	for (int i = 0; i < lights.size(); i++) {
		Eigen::Vector3f lightDirection = -(pointP - lights[i].first).normalized();
		sumColor += calcSingleColor(minimum_face, origin, lightDirection, lights[i].second, pointP);
	}
	for (int i = 0; i < dirLights.size(); i++) {
		sumColor += calcSingleColor(minimum_face, origin, get<0>(dirLights[i]), get<1>(dirLights[i]), pointP);
	}
	return Eigen::Vector3f(min(sumColor.x(), 1.f), min(sumColor.y(), 1.f), min(sumColor.z(), 1.f));
}
