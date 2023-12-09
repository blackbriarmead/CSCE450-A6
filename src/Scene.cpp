#include <iostream>

#include "GLSL.h"
#include "Scene.h"
#include "Particle.h"
#include "Cloth.h"
#include "Shape.h"
#include "Program.h"
#include "Terrain.h"
#include "Ocean.h"
#include "Boat.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0)
{
	this->t_start = getTime();

}

Scene::~Scene()
{
}

double Scene::getTime() {
	unsigned long milliseconds_since_epoch =
	std::chrono::system_clock::now().time_since_epoch() /
	std::chrono::milliseconds(1);

	double seconds = milliseconds_since_epoch / 1000.0;
	//std::cout << "seconds since epoch: " << seconds << std::endl;
	return(seconds);
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 1e-1; //1e-3 default time step normal range 1e-2 to 1e-4
	
	grav << 0.0, -9.8, 0.0;
	
	int rows = 100;
	int cols = 100;
	double mass = 1.0;
	double alpha = 0e-1; //0e-1 default
	double damping = 1e-5; //1e-5 default
	double pradius = 0.01; // Particle radius, used for collisions
	Vector3d x00(-0.25, 0.5, 0.0);
	Vector3d x01(0.25, 0.5, 0.0);
	Vector3d x10(-0.25, 0.5, -0.5);
	Vector3d x11(0.25, 0.5, -0.5);
	cloth = make_shared<Cloth>(rows, cols, x00, x01, x10, x11, mass, alpha, damping, pradius);
	
	sphereShape = make_shared<Shape>();
	sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	
	auto sphere = make_shared<Particle>(sphereShape);
	spheres.push_back(sphere);
	sphere->r = 0.1;
	sphere->x = Vector3d(0.0, 0.2, 0.0);
	const int terrain_size = 256;
	const int ocean_size = 128;
	terrain = make_shared<Terrain>(terrain_size, terrain_size,0.05,10.0,5.0);
	ocean = make_shared<Ocean>(ocean_size, ocean_size, 0.05f, 3.0, 0.5);
	boat = make_shared<Boat>(RESOURCE_DIR + "ship.obj", Vector3d(0.0,0.0,2.0));
}

void Scene::init()
{
	sphereShape->init();
	cloth->init();
	terrain->init();
	ocean->init();
	boat->init();
}

void Scene::tare()
{
	/*for(auto s : spheres) {
		s->tare();
	}
	cloth->tare();*/
}

void Scene::reset()
{
	/*t = 0.0;
	for(auto s : spheres) {
		s->reset();
	}
	cloth->reset();*/
}

void Scene::step()
{
	float prevT = t;
	t = 1.0*(getTime() - t_start);

	float dt = t - prevT;
	
	// Move the sphere
	if(!spheres.empty()) {
		auto s = spheres.front();
		s->x(2) = 0.5 * sin(0.5*t);
	}
	
	// Simulate the cloth
	//cloth->step(h, grav, spheres);
	ocean->step(t);
	boat->step(dt, ocean);
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple) const
{
	//glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	//glUniform3fv(programs.at(0)->getUniform("kdFront"), 1, Vector3f(0.0, 1.0, 0.0).data());
	//glUniform3fv(programs.at(0)->getUniform("kdBack"), 1, Vector3f(0.0, 1.0, 0.0).data());
	//for(auto s : spheres) {
		//s->draw(MV, prog);
	//}
	//cloth->draw(MV, prog);
	terrain->draw(MV, prog);
	ocean->draw(MV, prog);
	boat->draw(MV, prog, progSimple);
	
	GLSL::checkError(GET_FILE_LINE);
}
