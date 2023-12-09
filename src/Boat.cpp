#include "Boat.h"


#include <iostream>
#include <fstream>
#include <math.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

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
	this->rotation = Vector3d(0.0, 0.0, 0.0);
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
	//generate_fp(0.25, 0.05, 20, 20);
	generate_fp(25.0, 5.0, 20, 20);
}

const glm::mat4x4 Boat::getRotationMatrix() {
	Quaterniond rq = AngleAxisd(rotation[0]+3.1415926, Vector3d::UnitX())
	* AngleAxisd(rotation[1], Vector3d::UnitY())
	* AngleAxisd(rotation[2], Vector3d::UnitZ());
	
	glm::quat glm_quat = glm::quat(rq.coeffs()[0], rq.coeffs()[1], rq.coeffs()[2], rq.coeffs()[3]);

	return(glm::toMat4(glm_quat));
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
	const double BUOYANCY_CONSTANT = 2000.0;

	//go through every point
	std::vector<Vector3d> forces = std::vector<Vector3d>();
	std::vector<Vector3d> torques = std::vector<Vector3d>();
	float particle_mass = mass / fp_positions.size();
	for (int i = 0; i < fp_positions.size(); ++i) {
		Vector3d position = fp_positions.at(i) / 100.0;
		Vector3d normal = fp_normals.at(i);
		glm::vec4 glmpos = glm::vec4(position[0], position[1], position[2], 1);
		glm::vec4 glmnor = glm::vec4(normal[0], normal[1], normal[2], 0);
		glmpos = glmpos * this->getRotationMatrix();
		glmnor = glmnor * this->getRotationMatrix();
		position = Vector3d(glmpos.x, glmpos.y, glmpos.z)/glmpos.w;
		normal = Vector3d(glmnor.x, glmnor.y, glmnor.z);
		
		Vector3d worldPosition = this->position + position;
		Vector3d r = this->position - position; //vector from point to center of mass
		float waterHeight = ocean->getHeight(worldPosition[0], worldPosition[2]);

		float depth = waterHeight - worldPosition[1];

		Vector3d gravity = Vector3d(0.0, -9.8, 0.0);

		Vector3d grav_force = gravity * particle_mass;
		Vector3d buoyancy_force = Vector3d(0.0, 0.0, 0.0);
		//negative depth has no buoyancy force
		
		
		if (depth > 0) {
			//buoyancy_force = Vector3d(0.0, abs(depth * BUOYANCY_CONSTANT), 0.0) - 0.1 * this->velocity;
			buoyancy_force = normal * abs(particle_mass * depth * BUOYANCY_CONSTANT) - 10.0 * glm::pow(this->velocity.norm(),1.0) * this->velocity;
		}

		Vector3d total_force = grav_force + buoyancy_force;

		Vector3d torque = total_force.cross(r);

		if (i == 0) {
			//std::cout << "depth: " << depth << std::endl;
			//std::cout << "waterHeight: " << waterHeight << std::endl;
			//std::cout << "this->position: " << this->position << std::endl;

			

			//std::cout << "r: " << r << std::endl;
			
			//std::cout << "torque: " << torque << std::endl;

		}

		forces.push_back(total_force);
		torques.push_back(torque);

	}

	Vector3d summed_forces = Vector3d(0.0, 0.0, 0.0);

	for (Vector3d& force : forces) {
		summed_forces += force;
	}

	Vector3d summed_torques = Vector3d(0.0, 0.0, 0.0);

	for (Vector3d& torque : torques) {
		summed_torques += torque;
	}

	//std::cout << "summed torques: " << summed_torques << std::endl;

	summed_torques *= 0.01;

	const double DAMPING = 10.0;

	Vector3d damping_force = -this->velocity * DAMPING;

	summed_forces += damping_force;

	//summed_forces += Vector3d(5.0, 0.0, 0.0);

	summed_forces *= 0.01;

	this->velocity += dt * (summed_forces / mass);

	this->position += dt * velocity;

	this->angular_velocity += dt * summed_torques;
	this->angular_velocity *= 0.95;
	//Quaterniond angle_step = AngleAxisd(angular_velocity[0], Vector3d::UnitX())
	//	* AngleAxisd(angular_velocity[1], Vector3d::UnitY())
	//	* AngleAxisd(angular_velocity[2], Vector3d::UnitZ());
	this->rotation += dt * angular_velocity;
}

void Boat::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p, const shared_ptr<Program> progSimple)
{
	MV->pushMatrix();
	MV->multMatrix(getRotationMatrix());
	MV->translate(glm::vec3(position[0], position[1], position[2]));
	MV->scale(0.01f);
	glUniform3f(p->getUniform("kd"), 0.73f, 0.4137f, 0.0f);
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	boatMesh->draw(p);
	MV->popMatrix();
	MV->pushMatrix();
	MV->translate(glm::vec3(position[0], position[1], position[2]));
	
	if (false) {
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
	}
	

	MV->popMatrix();
}
