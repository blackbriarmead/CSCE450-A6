#include <iostream>
#include <fstream>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Cloth.h"
#include "Particle.h"
#include "Spring.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"

using namespace std;
using namespace Eigen;

Cloth::Cloth(int rows, int cols,
			 const Vector3d &x00,
			 const Vector3d &x01,
			 const Vector3d &x10,
			 const Vector3d &x11,
			 double mass,
			 double alpha,
			 double damping,
			 double pradius)
{
	assert(rows > 1);
	assert(cols > 1);
	assert(mass > 0.0);
	assert(alpha >= 0.0);
	assert(damping >= 0.0);
	assert(pradius >= 0.0);
	
	this->rows = rows;
	this->cols = cols;
	
	// Create particles
	int nVerts = rows*cols; // Total number of vertices
	double particleMass = mass / nVerts; //Spread mass of cloth evenly over all particles
	Vector3d u = x10 - x00;
	Vector3d du = u / cols;
	Vector3d v = x01 - x00;
	Vector3d dv = v / rows;
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			auto p = make_shared<Particle>();
			particles.push_back(p);
			p->r = pradius;
			p->m = particleMass;
			p->d = damping;
			
			p->x = x00 + i * dv + j * du;
			p->p = p->x;
			p->v = Eigen::Vector3d(0.0, 0.0, 0.0);
			
			p->fixed = (i == 0 && j == 0) || (i == rows - 1 && j == 0);//two corners should be fixed
		}
	}

	
	// Create x springs
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols-1; ++j) {
			springs.push_back(make_shared<Spring>(particles[i*cols+j],particles[i*cols+j+1],alpha));
		}
	}
	
	// Create y springs
	for (int i = 0; i < rows - 1; ++i) {
		for (int j = 0; j < cols; ++j) {
			springs.push_back(make_shared<Spring>(particles[i * cols + j], particles[(i+1) * cols + j], alpha));
		}
	}

	
	// Create shear springs
	for (int i = 0; i < rows - 1; ++i) {
		for (int j = 0; j < cols - 1; ++j) {
			springs.push_back(make_shared<Spring>(particles[i * cols + j], particles[(i+1) * cols + j + 1], alpha));
			springs.push_back(make_shared<Spring>(particles[i * cols + j + 1], particles[(i+1) * cols + j], alpha));
		}
	}
	
	// Create x bending springs
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols - 2; ++j) {
			springs.push_back(make_shared<Spring>(particles[i * cols + j], particles[i * cols + j + 2], alpha));
		}
	}

	// Create y bending springs
	for (int i = 0; i < rows - 2; ++i) {
		for (int j = 0; j < cols; ++j) {
			springs.push_back(make_shared<Spring>(particles[i * cols + j], particles[(i + 2) * cols + j], alpha));
		}
	}
	

	// Build vertex buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	posBuf.resize(nVerts*3);
	norBuf.resize(nVerts*3);
	updatePosNor();
	
	// Texture coordinates (don't change)
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			texBuf.push_back(i/(rows-1.0));
			texBuf.push_back(j/(cols-1.0));
		}
	}
	
	// Elements (don't change)
	for(int i = 0; i < rows-1; ++i) {
		for(int j = 0; j < cols; ++j) {
			int k0 = i*cols + j;
			int k1 = k0 + cols;
			// Triangle strip
			eleBuf.push_back(k0);
			eleBuf.push_back(k1);
		}
	}
}

Cloth::~Cloth()
{
}

void Cloth::tare()
{
	for(auto p : particles) {
		p->tare();
	}
}

void Cloth::reset()
{
	for(auto p : particles) {
		p->reset();
	}
	updatePosNor();
}

void Cloth::updatePosNor()
{
	// Position
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			int k = i*cols + j;
			Vector3d x = particles[k]->x; // updated position
			posBuf[3*k+0] = x(0);
			posBuf[3*k+1] = x(1);
			posBuf[3*k+2] = x(2);
		}
	}
	
	// Normal
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
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
			int k = i*cols + j;
			int ku0 = k - 1;
			int ku1 = k + 1;
			int kv0 = k - cols;
			int kv1 = k + cols;
			Vector3d x = particles[k]->x;
			Vector3d xu0, xu1, xv0, xv1, dx0, dx1, c;
			Vector3d nor(0.0, 0.0, 0.0);
			int count = 0;
			// Top-right triangle
			if(j != cols-1 && i != rows-1) {
				xu1 = particles[ku1]->x;
				xv1 = particles[kv1]->x;
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Top-left triangle
			if(j != 0 && i != rows-1) {
				xu1 = particles[kv1]->x;
				xv1 = particles[ku0]->x;
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Bottom-left triangle
			if(j != 0 && i != 0) {
				xu1 = particles[ku0]->x;
				xv1 = particles[kv0]->x;
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			// Bottom-right triangle
			if(j != cols-1 && i != 0) {
				xu1 = particles[kv0]->x;
				xv1 = particles[ku1]->x;
				dx0 = xu1 - x;
				dx1 = xv1 - x;
				c = dx0.cross(dx1);
				nor += c.normalized();
				++count;
			}
			nor /= count;
			nor.normalize();
			norBuf[3*k+0] = nor(0);
			norBuf[3*k+1] = nor(1);
			norBuf[3*k+2] = nor(2);
		}
	}

}

void Cloth::step(double h, const Vector3d &grav, const vector< shared_ptr<Particle> > spheres)
{
	//for each particles
	for (int i = 0; i < particles.size(); i++) {
		std::shared_ptr<Particle> particle = particles.at(i);
		if (particle->fixed) continue;

		//calculate the force on the particle
		auto fi = particle->m * grav - particle->d * particle->v;

		// I think h represents time step size?
		//calculate new velocity of particle from force and time step
		particle->v += (h / particle->m) * fi;
		//update previous pos before doing iteration
		particle->p = particle->x;
		//update particle position based on timestep and velocity
		particle->x += h * particle->v;
	}

	//for each constraint
	for (int i = 0; i < springs.size(); i++) {
		std::shared_ptr<Spring> spring = springs.at(i);
		Vector3d dx = spring->p1->x - spring->p0->x;
		double l = dx.norm();

		double C = l - spring->L;

		Vector3d grad_C0 = -(dx / l);
		Vector3d grad_C1 = -grad_C0;		

		//w is inverse of mass
		double w0 = 1.0 / spring->p0->m;
		double w1 = 1.0 / spring->p1->m;

		double lambdaj = -C / (w0 * (grad_C0.norm() * grad_C0.norm()) + w1 * (grad_C1.norm() * grad_C1.norm()) + spring->alpha / (h * h));
	
		if(!spring->p0->fixed) spring->p0->x += lambdaj * w0 * grad_C0;
		if(!spring->p1->fixed) spring->p1->x += lambdaj * w1 * grad_C1;
	
	}

	//do collision check
	for (int i = 0; i < particles.size(); i++) {
		std::shared_ptr<Particle> particle = particles.at(i);
		if (particle->fixed) continue;

		for (int j = 0; j < spheres.size(); j++) {
			std::shared_ptr<Particle> sphere = spheres.at(j);
			//Vector from sphere to particle
			Vector3d x = particle->x - sphere->x;
			double distance = x.norm();
			//detect collision
			double margin = distance - sphere->r - particle->r;
			if (margin < 0) {
				//move margin along line away from sphere
				Vector3d toMove = x/x.norm() * -margin;
				particle->x += toMove;
			}
		}
	}

	//for each particle
	for (int i = 0; i < particles.size(); i++) {
		std::shared_ptr<Particle> particle = particles.at(i);
		if (particle->fixed) continue;

		//set particle's velocity to return to previous point
		particle->v = (1/h) * (particle->x - particle->p);
	}

	// Update position and normal buffers
	updatePosNor();
}

void Cloth::init()
{
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	
	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	
	glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_STATIC_DRAW);
	
	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size()*sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	
	assert(glGetError() == GL_NO_ERROR);
}

void Cloth::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
	// Draw mesh
	glUniform3f(p->getUniform("kd"), 0.894f, 0.882f, 0.792f);
	//glUniform3f(p->getUniform("kdBack"),  0.776f, 0.843f, 0.835f);
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	int h_tex = p->getAttribute("aTex");
	if(h_tex >= 0) {
		glEnableVertexAttribArray(h_tex);
		glBindBuffer(GL_ARRAY_BUFFER, texBufID);
		glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	for(int i = 0; i < rows; ++i) {
		glDrawElements(GL_TRIANGLE_STRIP, 2*cols, GL_UNSIGNED_INT, (const void *)(2*cols*i*sizeof(unsigned int)));
	}
	if(h_tex >= 0) {
		glDisableVertexAttribArray(h_tex);
	}
	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}
