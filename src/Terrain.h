#pragma once
#ifndef Terrain_H
#define Terrain_H

#include <map>
#include <string>
#include <Eigen/Dense>

#define GLEW_STATIC
#include <GL/glew.h>
#include <vector>
#include <memory>


class MatrixStack;
class Program;

class Terrain
{
public:
	Terrain(double resolution, int width, int height);
	virtual ~Terrain();
	void init();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void updatePosNor();

private:
	int width;
	int height;
	double resolution;
	Eigen::Vector3d centerOffset;

	std::vector<float> heightMap;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

#endif
