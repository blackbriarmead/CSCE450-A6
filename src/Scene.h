#pragma once
#ifndef Scene_H
#define Scene_H

#include <vector>
#include <memory>
#include <string>
#include <chrono>
#include <iostream>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Cloth;
class Particle;
class MatrixStack;
class Program;
class Shape;
class Terrain;
class Ocean;
class Boat;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	double getTime();
	void init();
	void tare();
	void reset();
	void step();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple) const;
	
	double getTime() const { return t; }
	std::vector<std::shared_ptr<Program>> programs;
	std::shared_ptr<Boat> boat;
	
private:
	double t;
	double t_start;
	double h;
	Eigen::Vector3d grav;
	
	std::shared_ptr<Shape> sphereShape;
	std::shared_ptr<Cloth> cloth;
	std::vector< std::shared_ptr<Particle> > spheres;
	std::shared_ptr<Terrain> terrain;
	std::shared_ptr<Ocean> ocean;
	
	
};

#endif
