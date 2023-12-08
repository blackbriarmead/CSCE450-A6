#pragma once
#ifndef Boat_H
#define Boat_H

#include <map>
#include <string>
#include <Eigen/Dense>

#define GLEW_STATIC
#include <GL/glew.h>
#include <vector>
#include <memory>
#include <iostream>


class MatrixStack;
class Program;
class Shape;
class Ocean;

class Boat
{
public:
	Boat(std::string meshName, Eigen::Vector3d startPos);
	virtual ~Boat();
	void init();
	void step(double dt, std::shared_ptr<Ocean> ocean);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p, const std::shared_ptr<Program> progSimple) const;
	void generate_fp(float length, float radius, float lengthwise_points, float horizontal_points);

	Eigen::Vector3d position;
	Eigen::Vector3d velocity;
	Eigen::Vector3d angular_velocity;

	//fp stands for flotation particle
	//these particles interact with the ocean mesh to calculate the force they experience
	//this force is summed up for all of these particles, then applied to the boat as a whole
	std::vector<Eigen::Vector3d> fp_positions;
	std::vector<Eigen::Vector3d> fp_normals;
	float mass;

private:
	std::shared_ptr<Shape> boatMesh;
};

#endif
