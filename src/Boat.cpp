#include "Boat.h"


#include <iostream>
#include <fstream>
#include <math.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Particle.h"
#include "Spring.h"
#include "Shape.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include "PerlinNoise.hpp"
#include "Ocean.h"

using namespace std;
using namespace Eigen;

Boat::Boat(std::string meshName, Eigen::Vector3d startPos)
{
	boatMesh = make_shared<Shape>();
	boatMesh->loadMesh(meshName);
	this->position = startPos;
	this->velocity = Vector3d(0.0, 0.0, 0.0);
	fp_positions = std::vector<Vector3d>();
	fp_normals = std::vector<Vector3d>();
	mass = 1.0;
}

Boat::~Boat()
{
}

void Boat::init()
{
	boatMesh->init();
	generate_fp(25.0, 5.0, 20, 20);
}


//Generates the points that are used to calculate the forces applied to the boat
void Boat::generate_fp(float length, float radius, float lengthwise_points, float horizontal_points)
{
	for (int zi = 0; zi < lengthwise_points; ++zi) {
		float z = -0.5 * length + length * (zi / (lengthwise_points - 1));

		for (int ti = 0; ti < horizontal_points; ++ti) {
			float theta = 3.1415926 + 3.1415926 * (ti / (horizontal_points - 1));

			Vector3d fp_pos = Vector3d(z, radius * sin(theta), radius * cos(theta));
			Vector3d fp_nor = Vector3d(0, -fp_pos[1], -fp_pos[2]);
			fp_nor.normalize();

			fp_positions.push_back(fp_pos);
			fp_normals.push_back(fp_nor);
		}
	}
}

void Boat::step(double dt, shared_ptr<Ocean> ocean) {
	//go through every point
	std::vector<Vector3d> forces = std::vector<Vector3d>();
	float particle_mass = mass / fp_positions.size();
	for (int i = 0; i < fp_positions.size(); ++i) {
		Vector3d position = fp_positions.at(i);
		Vector3d normal = fp_normals.at(i);
		
		Vector3d worldPosition = this->position + position;
		float waterHeight = ocean->getHeight(worldPosition[0], worldPosition[2]);

		float depth = waterHeight - worldPosition[1];

		Vector3d gravity = Vector3d(0.0, -9.8, 0.0);

		Vector3d grav_force = gravity * particle_mass;
		Vector3d buoyancy_force = Vector3d(0.0, 0.0, 0.0);
		//negative depth has no buoyancy force
		if (depth > 0) {

		}

	}


	position += dt * velocity;
}

void Boat::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p, const shared_ptr<Program> progSimple) const
{
	MV->pushMatrix();
	MV->translate(glm::vec3(position[0], position[1], position[2]));
	MV->scale(0.01f);
	glUniform3f(p->getUniform("kd"), 0.73f, 0.4137f, 0.0f);
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	boatMesh->draw(p);
	glPointSize(10.0f);
	glBegin(GL_POINTS);
	glColor3f(1.0, 0.5, 0.0);
	for (const Vector3d& position : fp_positions) {
		glVertex3f(position[0], position[1], position[2]);
	}
	glEnd();

	glBegin(GL_LINES);
	glColor3f(0.0, 1.0, 0.0);
	for (int i = 0; i < fp_positions.size(); ++i) {
		Vector3d position = fp_positions.at(i);
		Vector3d normal = fp_normals.at(i);
		glVertex3f(position[0], position[1], position[2]);
		glVertex3f(position[0] + normal[0], position[1] + normal[1], position[2] + normal[2]);
	}
	glEnd();

	MV->popMatrix();
}
