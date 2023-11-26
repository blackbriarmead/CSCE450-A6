#include "Terrain.h"


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

Terrain::Terrain(double resolution,int width, int height)
{
	assert(width > 1);
	assert(height > 1);
	assert(resolution > 0.0);

	this->width = width;
	this->height = height;
	this->resolution = resolution;
	this->centerOffset = Vector3d(width * resolution / 2, 1, height * resolution / 2);

	int nVerts = width * height;

	//Generate Terrain using noise generator
	siv::PerlinNoise perlin;
	for (int row = 0; row < height; ++row) {
		for (int col = 0; col < width; ++col) {
			heightMap.push_back(perlin.normalizedOctave2D(row * resolution, col * resolution, 10));
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

Terrain::~Terrain()
{
}

void Terrain::updatePosNor()
{
	// Position
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			int k = i * width + j;
			posBuf[3 * k + 0] = j * resolution - centerOffset[0];
			//posBuf[3 * k + 1] = -1.0f;
			posBuf[3 * k + 1] = 1*heightMap.at(k) - centerOffset[1];
			posBuf[3 * k + 2] = i * resolution - centerOffset[2];
		}
	}

	std::cout << posBuf.size() << std::endl;
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
			Vector3d x(j / resolution, heightMap[k], i / resolution);
			Vector3d xu0, xu1, xv0, xv1, dx0, dx1, c;
			Vector3d nor(0.0, 0.0, 0.0);
			int count = 0;
			// Top-right triangle
			if (j != width - 1 && i != height - 1) {
				xu1 = Vector3d((j+1) / resolution, heightMap[k+1], i / resolution);
				xv1 = Vector3d((j + 0) / resolution, heightMap[k+width], (i+1) / resolution);
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Top-left triangle
			if (j != 0 && i != height - 1) {
				xu1 = Vector3d((j + 0) / resolution, heightMap[k + width], (i + 1) / resolution);
				xv1 = Vector3d((j - 1) / resolution, heightMap[k-1], i / resolution);
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Bottom-left triangle
			if (j != 0 && i != 0) {
				xu1 = Vector3d((j - 1) / resolution, heightMap[k - 1], i / resolution);
				xv1 = Vector3d((j + 0) / resolution, heightMap[k - width], (i - 1) / resolution);
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Bottom-right triangle
			if (j != width - 1 && i != 0) {
				xu1 = Vector3d((j + 0) / resolution, heightMap[k - width], (i - 1) / resolution);
				xv1 = Vector3d((j + 1) / resolution, heightMap[k + 1], i / resolution);
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			nor /= count;
			nor.normalize();
			norBuf[3 * k + 0] = nor(0);
			norBuf[3 * k + 1] = nor(1);
			norBuf[3 * k + 2] = nor(2);
		}
	}

}

void Terrain::init()
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

void Terrain::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
	// Draw mesh
	glUniform3f(p->getUniform("kdFront"), 0.894f, 0.882f, 0.792f);
	glUniform3f(p->getUniform("kdBack"), 0.776f, 0.843f, 0.835f);
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);
	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);
	int h_tex = p->getAttribute("aTex");
	if (h_tex >= 0) {
		glEnableVertexAttribArray(h_tex);
		glBindBuffer(GL_ARRAY_BUFFER, texBufID);
		glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, (const void*)0);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	for (int i = 0; i < height; ++i) {
		glDrawElements(GL_TRIANGLE_STRIP, 2 * width, GL_UNSIGNED_INT, (const void*)(2 * width * i * sizeof(unsigned int)));
	}
	if (h_tex >= 0) {
		glDisableVertexAttribArray(h_tex);
	}
	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}
