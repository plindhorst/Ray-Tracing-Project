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

void Flyscene::simulate(GLFWwindow *window) {
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

void Flyscene::createDebugRay(const Eigen::Vector2f &mouse_pos) {
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

  // create spherical lights out of point lights
  vector<Eigen::Vector3f> lightsDup = lights;
  for (int i = 0; i < lightsDup.size(); i++) {
	  sphericalLight(lightsDup[i], 0.15, 15);
  }

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


Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f &origin,
                                   Eigen::Vector3f &dest) {
  // just some fake random color per pixel until you implement your ray tracing
  // remember to return your RGB values as floats in the range [0, 1]!!!
  return Eigen::Vector3f(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX,
                         rand() / (float)RAND_MAX);
}

float Flyscene::intersection(Eigen::Vector3f& origin, Eigen::Vector3f& dest, Tucano::Face& face) {
	return -1;
}

//Call function in calculate color, if it returns true => pixel should be black/ambiant. If it returns false => the pixel should have a color.
bool Flyscene::shadow(Eigen::Vector3f& loc, Eigen::Vector3f& lightLoc) {
	Eigen::Vector3f lightDirection = loc - lightLoc;
	Tucano::Face current_face;
	for (int i = 0; i < Flyscene::mesh.getNumberOfFaces(); i++) {
		current_face = mesh.getFace(i);
		if (intersection(loc, lightDirection, current_face) != -1) {
			return true;
		}
	}
	return false;
}

//Call function in rayTraceScene before pixel loop
void Flyscene::sphericalLight(Eigen::Vector3f& lightLoc, float radius, int nLightpoints) {
	for (int i = 0; i < nLightpoints; i++) {
		float x = (-radius + (rand() / (RAND_MAX / (radius * 2))));
		float y = (-radius + (rand() / (RAND_MAX / (radius * 2))));
		float z = (-radius + (rand() / (RAND_MAX / (radius * 2))));
		lights.push_back(Eigen::Vector3f(lightLoc.x() + x, lightLoc.y() + y, lightLoc.z() + z));
	}
}


//Call function in calcColor
Eigen::Vector3f Flyscene::calcSingleColor(Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f& lightLoc, Eigen::Vector3f& pointP) {
	
	if (shadow(pointP, lightLoc)) {
		return Eigen::Vector3f(0.0, 0.0, 0.0);
	}

	Tucano::Material::Mtl material = materials[minimum_face.material_id];
	Eigen::Vector3f ka = material.getAmbient();
	Eigen::Vector3f kd = material.getDiffuse();
	Eigen::Vector3f ks = material.getSpecular();
	float shininess = material.getShininess();
	float refraction_index = material.getOpticalDensity();
	float transparency = material.getDissolveFactor();

	Eigen::Vector3f normal = minimum_face.normal.normalized();
	Eigen::Vector3f lightDirection = (pointP - lightLoc).normalized();
	Eigen::Vector3f lightReflection = (lightDirection - 2 * (normal.dot(lightDirection)) * normal).normalized();
	Eigen::Vector3f eyeDirection = (origin - pointP).normalized();
	Eigen::Vector3f light_intensity = Eigen::Vector3f(1.0, 1.0, 1.0);
	
	Eigen::Vector3f ambient = Eigen::Vector3f(light_intensity.x() * ka.x(), light_intensity.y() * ka.y(), light_intensity.z() * ka.z());
	Eigen::Vector3f diffuse = Eigen::Vector3f(light_intensity.x() * ks.x(), light_intensity.y() * ks.y(), light_intensity.z() * ks.z()) * std::max(lightDirection.dot(normal), 0.f);
	Eigen::Vector3f specular = Eigen::Vector3f(light_intensity.x() * kd.x(), light_intensity.y() * kd.y(), light_intensity.z() * kd.z()) * std::max(std::pow(lightReflection.dot(eyeDirection), shininess), 0.f);

	Eigen::Vector3f color = ambient + diffuse + specular;
	return color;
}

//Call function in traceRay
Eigen::Vector3f Flyscene::calcColor(Tucano::Face minimum_face, Eigen::Vector3f& origin, Eigen::Vector3f& pointP) {
	Eigen::Vector3f sumColor = Eigen::Vector3f(0.0, 0.0, 0.0);
	for (int i = 0; i < lights.size(); i++) {
		sumColor = sumColor + calcSingleColor(minimum_face, origin, lights[i], pointP);
	}
	return sumColor / lights.size();
}
