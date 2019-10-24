#include <algorithm>

#include "flyscene.hpp"
#include <GLFW/glfw3.h>

void Flyscene::initialize(int width, int height) {
	// initiliaze the Phong Shading effect for the Opengl Previewer
	phong.initialize();

	// set the camera's projection matrix
	flycamera.setPerspectiveMatrix(60.0, width / (float)height, 0.1f, 100.0f);
	flycamera.setViewport(Eigen::Vector2f((float)width, (float)height));

	// load the OBJ file and materials
	Tucano::MeshImporter::loadObjFile(mesh, materials,
		"resources/models/cube.obj");


	// normalize the model (scale to unit cube and center at origin)
	mesh.normalizeModelMatrix();

	// pass all the materials to the Phong Shader
	for (int i = 0; i < materials.size(); ++i)
		phong.addMaterial(materials[i]);



	// set the color and size of the sphere to represent the light sources
	// same sphere is used for all sources
	lightrep.setColor(Eigen::Vector4f(1.0, 1.0, 0.0, 1.0));
	lightrep.setSize(0.15);

	// create a first ray-tracing light source at some random position
	lights.push_back(Eigen::Vector3f(-1.0, 1.0, 1.0));

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
	scene_light.viewMatrix()->translate(-lights.back());

	// render the scene using OpenGL and one light source
	phong.render(mesh, flycamera, scene_light);

	// render the ray and camera representation for ray debug
	ray.render(flycamera, scene_light);
	camerarep.render(flycamera, scene_light);

	// render ray tracing light sources as yellow spheres
	for (int i = 0; i < lights.size(); ++i) {
		lightrep.resetModelMatrix();
		lightrep.modelMatrix()->translate(lights[i]);
		lightrep.render(flycamera, scene_light);
	}

	// render coordinate system at lower right corner
	flycamera.renderAtCorner();
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

	Eigen::Vector3f origin = flycamera.getCenter();
	Eigen::Vector3f debugColor = traceRay(origin, dir);
	std::cout << "Debug Ray Color: " << debugColor.transpose() << std::endl;
}

void Flyscene::raytraceScene(int width, int height) {
	std::cout << "ray tracing ..." << std::endl;

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

	// origin of the ray is always the camera center
	Eigen::Vector3f origin = flycamera.getCenter();
	Eigen::Vector3f screen_coords;

	// for every pixel shoot a ray from the origin through the pixel coords
	for (int j = 0; j < image_size[1]; ++j) {
		for (int i = 0; i < image_size[0]; ++i) {
			// create a ray from the camera passing through the pixel (i,j)
			screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
			// launch raytracing for the given ray and write result to pixel data
			pixel_data[i][j] = traceRay(origin, screen_coords);
		}
	}

	// write the ray tracing result to a PPM image
	Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
	std::cout << "ray tracing done! " << std::endl;
}


Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f& origin, Eigen::Vector3f& dest) {
	// Parameters to keep track of current faces and the closest face
	float minimum_distance = FLT_MAX;
	Tucano::Face minimum_face;
	float current_distance;
	Tucano::Face current_face;

	// Loop through all faces, calculate the distance and update the minimum_face whenever necessary
	for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
		current_face = mesh.getFace(i);
		current_distance = calculateDistance(origin, dest, current_face);
		if (0 <= current_distance && current_distance < minimum_distance) {
			minimum_distance = current_distance;
			minimum_face = current_face;
		}
	}
	Eigen::Vector3f color;

	// Test if the ray intersected with a face, if so: calculate the color
	if (minimum_distance != FLT_MAX) {
		color = calculateColor(minimum_distance, minimum_face, origin, dest);
	}
	else {
		color = Eigen::Vector3f(0.9, 0.9, 0.9);
	}
	return color;
}

float Flyscene::calculateDistance(Eigen::Vector3f& origin, Eigen::Vector3f& dest, Tucano::Face face) {
	//get vertices and normals of the face
	Eigen::Affine3f modelMatrix = mesh.getShapeModelMatrix();

	Eigen::Vector4f vec0 = mesh.getVertex(face.vertex_ids[0]);
	Eigen::Vector4f vec1 = mesh.getVertex(face.vertex_ids[1]);
	Eigen::Vector4f vec2 = mesh.getVertex(face.vertex_ids[2]);

	Eigen::Vector3f vert0 = modelMatrix * (vec0.head<3>()/vec0.w());
	Eigen::Vector3f vert1 = modelMatrix * (vec1.head<3>()/vec1.w());
	Eigen::Vector3f vert2 = modelMatrix * (vec2.head<3>()/vec2.w());



	//get Normal of the face
	Eigen::Vector3f facenormal = face.normal;

	//Return false if triangle and direction of ray are the same
	if (facenormal.dot(dest) == 0) {
		return (float)-1;
	}
	//Calculate the distance between the plane and the origin (not the camera)
	float distancePlane = facenormal.dot(vert0.normalized());
	//Ray = origin + t*distance

	float orthProjectionDest = distancePlane - origin.dot(facenormal);
	float t = orthProjectionDest / (dest.dot(facenormal));
	Eigen::Vector3f PointP = origin + t * dest;

	//Inside-out test
	Eigen::Vector3f edge0 = vert1 - vert0;
	Eigen::Vector3f edge1 = vert2 - vert1;
	Eigen::Vector3f edge2 = vert0 - vert2;

	Eigen::Vector3f Inner0 = PointP - vert0;
	Eigen::Vector3f Inner1 = PointP - vert1;
	Eigen::Vector3f Inner2 = PointP - vert2;

	float area0 = facenormal.dot(edge0.cross(Inner0));
	float area1 = facenormal.dot(edge1.cross(Inner1));
	float area2 = facenormal.dot(edge2.cross(Inner2));

	//If any area is smaller or equal to zero, point is on the wrong side of the edge. 
	if (area0 < 0 || area1 < 0 || area2 < 0) {
		return (float)-1;
	}
	else {
		Eigen::Vector3f distVector = t * dest;
		float distance = ::sqrtf(distVector.x()*distVector.x() + distVector.y()*distVector.y() + distVector.z()*distVector.z());
		return distance;
	}
}

Eigen::Vector3f Flyscene::calculateColor(float minimum_distance, Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f& dest) {
	Tucano::Material::Mtl material = materials[minimum_face.material_id];
	Eigen::Vector3f ka = material.getAmbient();
	Eigen::Vector3f kd = material.getDiffuse();
	Eigen::Vector3f ks = material.getSpecular();
	float shininess = material.getShininess();
	float refraction_index = material.getOpticalDensity();
	float transparency = material.getDissolveFactor();

	Eigen::Affine3f viewMatrix = flycamera.getViewMatrix();
	Eigen::Matrix4f projectionMatrix = flycamera.getProjectionMatrix();
	Eigen::Affine3f modelMatrix = mesh.getShapeModelMatrix();
	Eigen::Affine3f lightViewMatrix = scene_light.getViewMatrix();

	Eigen::Vector3f normal = minimum_face.normal.normalized();
	Eigen::Vector3f lightDirection = (viewMatrix * lightViewMatrix.inverse() * Eigen::Vector3f(0.0, 0.0, 1.0)).normalized();
	Eigen::Vector3f lightReflection = (-lightDirection - 2 * (normal.dot(-lightDirection)) * normal).normalized();
	Eigen::Vector3f eyeDirection = (origin - dest).normalized();

	Eigen::Vector3f light_intensity = Eigen::Vector3f(1.0, 1.0, 1.0);

	Eigen::Vector3f ambient = Eigen::Vector3f(light_intensity.x() * ka.x(), light_intensity.y() * ka.y(), light_intensity.z() * ka.z());
	Eigen::Vector3f diffuse = Eigen::Vector3f(light_intensity.x() * ks.x(), light_intensity.y() * ks.y(), light_intensity.z() * ks.z()) * std::max(lightDirection.dot(normal), 0.f);
	Eigen::Vector3f specular = Eigen::Vector3f(light_intensity.x() * kd.x(), light_intensity.y() * kd.y(), light_intensity.z() * kd.z()) * std::max(std::pow(lightReflection.dot(eyeDirection), shininess), 0.f);

	Eigen::Vector3f color = ambient + diffuse + specular;
	return color;
}
