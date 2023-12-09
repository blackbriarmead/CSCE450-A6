#include "Ocean.h"


#include <iostream>
#include <fstream>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Particle.h"
#include "Spring.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include "PerlinNoise.hpp"

using namespace std;
using namespace Eigen;

Ocean::Ocean(int width, int height, double resolution, double scale, double heightMul)
{
	assert(width > 1);
	assert(height > 1);
	assert(resolution > 0.0);

	this->width = width;
	this->height = height;
	this->resolution = resolution;
	this->scale = scale;
	this->heightMul = heightMul;
	this->centerOffset = Vector3d(width * resolution / 2, 0.1*heightMul, height * resolution / 2);

	int nVerts = width * height;

	//Generate Ocean using noise generator
	siv::PerlinNoise perlin;
	perlin.reseed(0);
	for (int row = 0; row < height; ++row) {
		for (int col = 0; col < width; ++col) {
			heightMap.push_back(perlin.normalizedOctave3D(row * resolution/scale, col * resolution/scale, 0, 10,0.5));
		}
	}

	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	posBuf.resize(nVerts * 3);
	norBuf.resize(nVerts * 3);
	updatePosNor();

	// Texture coordinates (don't change)
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			texBuf.push_back(i / (height - 1.0));
			texBuf.push_back(j / (width - 1.0));
		}
	}

	// Elements (don't change)
	for (int i = 0; i < height - 1; ++i) {
		for (int j = 0; j < width; ++j) {
			int k0 = i * width + j;
			int k1 = k0 + width;
			// Triangle strip
			eleBuf.push_back(k0);
			eleBuf.push_back(k1);
		}
	}
}

Ocean::~Ocean()
{
}

void Ocean::updatePosNor()
{
	// Position
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			int k = i * width + j;
			posBuf[3 * k + 0] = j * resolution - centerOffset[0];
			//posBuf[3 * k + 1] = -1.0f;
			posBuf[3 * k + 1] = this->heightMul*heightMap.at(k) - centerOffset[1];
			posBuf[3 * k + 2] = i * resolution - centerOffset[2];
		}
	}

	//std::cout << posBuf.size() << std::endl;
	/*for (float f : heightMap) {
		std::cout << f << std::endl;
	}*/

	// Normal
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			// Each particle has four neighbors
			//
			//      v1
			//     / | \
			// u0 /__|__\ u1
			//    \  |  /
			//     \ | /
			//      v0
			//
			// Use these four triangles to compute the normal
			int k = i * width + j;
			int ku0 = k - 1;
			int ku1 = k + 1;
			int kv0 = k - width;
			int kv1 = k + width;
			Vector3d x(posBuf[3*k], posBuf[3 * k + 1], posBuf[3 * k + 2]);
			Vector3d xu0, xu1, xv0, xv1, dx0, dx1, c;
			Vector3d nor(0.0, 0.0, 0.0);
			int count = 0;
			// Top-right triangle
			if (j != width - 1 && i != height - 1) {
				xu1 = Vector3d(posBuf[3 * ku1], posBuf[3 * ku1 + 1], posBuf[3 * ku1 + 2]);
				xv1 = Vector3d(posBuf[3 * kv1], posBuf[3 * kv1 + 1], posBuf[3 * kv1 + 2]);
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Top-left triangle
			if (j != 0 && i != height - 1) {
				xu1 = Vector3d(posBuf[3 * kv1], posBuf[3 * kv1 + 1], posBuf[3 * kv1 + 2]);
				xv1 = Vector3d(posBuf[3 * ku0], posBuf[3 * ku0 + 1], posBuf[3 * ku0 + 2]);
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Bottom-left triangle
			if (j != 0 && i != 0) {
				xu1 = Vector3d(posBuf[3 * ku0], posBuf[3 * ku0 + 1], posBuf[3 * ku0 + 2]);
				xv1 = Vector3d(posBuf[3 * kv0], posBuf[3 * kv0 + 1], posBuf[3 * kv0 + 2]);
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Bottom-right triangle
			if (j != width - 1 && i != 0) {
				xu1 = Vector3d(posBuf[3 * kv0], posBuf[3 * kv0 + 1], posBuf[3 * kv0 + 2]);
				xv1 = Vector3d(posBuf[3 * ku1], posBuf[3 * ku1 + 1], posBuf[3 * ku1 + 2]);
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			nor /= count;
			nor.normalize();
			//nor = -nor;
			norBuf[3 * k + 0] = nor(0);
			norBuf[3 * k + 1] = nor(1);
			norBuf[3 * k + 2] = nor(2);
		}
	}

}

void Ocean::init()
{
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size() * sizeof(float), &texBuf[0], GL_STATIC_DRAW);

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	assert(glGetError() == GL_NO_ERROR);
}

void Ocean::step(double t) {
	//Generate Ocean using noise generator
	siv::PerlinNoise perlin;
	perlin.reseed(0);
	heightMap.clear();
	for (int row = 0; row < height; ++row) {
		for (int col = 0; col < width; ++col) {
			heightMap.push_back(perlin.normalizedOctave3D(row * resolution / scale, col * resolution / scale, t/20.0, 10, 0.5));
		}
	}

	updatePosNor();
}

//Given x and z in world-space coordinates, return y, which is the height of the ocean at that point
float Ocean::getHeight(double x, double z) {

	// x = j * resolution - centerOffset[0]
	// centerOffset[0] + x = j * resolution
	// j = (centerOffset[0] + x)/resolution
	// z - i * resolution - centerOffset[2]
	// i = (centerOffset[2] + z)/resolution
	double i = (centerOffset[2] + z) / resolution;
	double j = (centerOffset[0] + x) / resolution;	

	//now that we have the coordinates, maybe just start by taking floor. better solution would be 2d lerp.
	int row = floor(i);
	int col = floor(j);

	row = max(0, min(height - 1, row));
	col = max(0, min(width - 1, col));
	//world = this->heightMul* heightMap.at(k) - centerOffset[1]
	return(this->heightMul * heightMap.at(row * width + col) - centerOffset[1]);
}

void Ocean::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
	// Draw mesh
	glUniform3f(p->getUniform("kd"), 0.0f, 0.0f, 1.0f);;
	//glUniform3f(p->getUniform("kdBack"), 0.776f, 0.843f, 0.835f);
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	GLSL::checkError(GET_FILE_LINE);
	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);
	GLSL::checkError(GET_FILE_LINE);
	int h_nor = p->getAttribute("aNor");
	GLSL::checkError(GET_FILE_LINE);
	glEnableVertexAttribArray(h_nor);
	GLSL::checkError(GET_FILE_LINE);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	GLSL::checkError(GET_FILE_LINE);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	GLSL::checkError(GET_FILE_LINE);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);
	GLSL::checkError(GET_FILE_LINE);
	int h_tex = p->getAttribute("aTex");
	if (h_tex >= 0) {
		glEnableVertexAttribArray(h_tex);
		glBindBuffer(GL_ARRAY_BUFFER, texBufID);
		glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, (const void*)0);
	}
	GLSL::checkError(GET_FILE_LINE);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	for (int i = 0; i < height; ++i) {
		glDrawElements(GL_TRIANGLE_STRIP, 2 * width, GL_UNSIGNED_INT, (const void*)(2 * width * i * sizeof(unsigned int)));
	}
	if (h_tex >= 0) {
		glDisableVertexAttribArray(h_tex);
	}
	GLSL::checkError(GET_FILE_LINE);
	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
	GLSL::checkError(GET_FILE_LINE);
}
