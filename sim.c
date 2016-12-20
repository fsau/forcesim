// Copyright (c) 2016, Franco Sauvisky
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//    derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define dim 2
#define maxParticles 10
#define maxForces 5
#define maxModifiers 5

#define maxTimeStep 1e-2
#define minTimeStep 1e-6
#define dt0 1e-4

#include "sim_internal.c"

#define tmax 10

void grav_force(double *, struct Particle, struct Particle);
void eletr_force(double *, struct Particle, struct Particle);
void mag_ext_force(double *, struct Particle, struct Particle);
void eletr_ext_force(double *, struct Particle, struct Particle);
void drag_force(double *, struct Particle, struct Particle);
void square_wall(struct Particle *, struct Particle *);
void sphere_collision(struct Particle *, struct Particle *);

int
main(void)
{
	struct point_system psystem = {
		.Bodies = {
			{
				.pos = {0,0},
				.vel = {0,0},
				.mass = 10,
				.radius = 0.5,
				.charge = 10,
				.fixed = false
			},
			{
				.pos = {3.9,4.1},
				.vel = {-8,-2},
				.mass = 10,
				.radius = 0.1,
				.charge = -10,
				.fixed = false
			},
			{
				.pos = {3,-3}, 
				.vel = {0,4},
				.mass = 10,
				.radius = 0.3,
				.charge = -10,
				.fixed = false
			},
			{
				.pos = {-3,3}, 
				.vel = {0,-4},
				.mass = 20,
				.radius = 0.3,
				.charge = -0,
				.fixed = false
			},
		},
		.Forces = {
			grav_force,
			eletr_force,
			mag_ext_force,
			eletr_ext_force,
			drag_force,
		},
		.Modifiers = {
			square_wall,
			sphere_collision,
		},
	};

	config_system(&psystem);

	printf("%d ", psystem.particleN);
	print_system(psystem);

	double dt = dynamic_dt(psystem);
	for(; psystem.time < tmax;)
	{
		dt = dynamic_dt(psystem);
		calc_system(&psystem, dt);
		print_system(psystem);
	}
}

#define grav_c 1
void
grav_force(double *force, struct Particle particle_0, struct Particle particle_1)
{
	if(particle_0.id == particle_1.id || particle_0.fixed) return;

	double r_01[dim];
	
	sum_vec(r_01, particle_1.pos, particle_0.pos, -1);

	double dist = r_abs(r_01);
	
	for(int i = 0; i < dim; i++)
	{
		force[i] += grav_c * particle_0.mass * particle_1.mass * r_01[i] / pow(dist, 3);
	}
}

#define eletr_c 2
void
eletr_force(double *force, struct Particle particle_0, struct Particle particle_1)
{
	if(particle_0.id == particle_1.id || particle_0.fixed) return;

	double r_01[dim];
	
	sum_vec(r_01, particle_1.pos, particle_0.pos, -1);

	double dist = r_abs(r_01);
	
	for(int i = 0; i < dim; i++)
	{
		force[i] -= eletr_c * particle_0.charge * particle_1.charge * r_01[i] / pow(dist, 3);
	}
}

#define mag_ext_c 1
void
mag_ext_force(double *force, struct Particle particle_0, struct Particle particle_1)
{
	if(particle_0.id != particle_1.id || particle_0.fixed) return;

	force[0] += mag_ext_c * particle_1.charge * particle_0.vel[1];
	force[1] -= mag_ext_c * particle_1.charge * particle_0.vel[0];
}

#define eletr_ext_c 20
void
eletr_ext_force(double *force, struct Particle particle_0, struct Particle particle_1)
{
	if(particle_0.id != particle_1.id || particle_0.fixed) return;

	force[0] += eletr_ext_c * particle_1.charge;
}

#define drag_c 2
void
drag_force(double *force, struct Particle particle_0, struct Particle particle_1)
{
	if(particle_0.id != particle_1.id || particle_0.fixed) return;

	for(int i = 0; i < dim; i++)
	{
		force[i] -= drag_c*particle_0.vel[i];
	}
}

#define topw 5
#define bottomw -5
#define rightw 5
#define leftw -5
void 
square_wall(struct Particle *particle_0, struct Particle *particle_1)
{
	if((*particle_0).id != (*particle_1).id) return;

	if((*particle_0).pos[1] + (*particle_0).radius > topw && (*particle_0).vel[1] > 0)
		(*particle_0).vel[1] *= -1;

	if((*particle_0).pos[0] + (*particle_0).radius > rightw && (*particle_0).vel[0] > 0)
		(*particle_0).vel[0] *= -1;

	if((*particle_0).pos[1] - (*particle_0).radius < bottomw && (*particle_0).vel[1] < 0)
		(*particle_0).vel[1] *= -1;

	if((*particle_0).pos[0] - (*particle_0).radius < leftw && (*particle_0).vel[0] < 0)
		(*particle_0).vel[0] *= -1;
}

void
sphere_collision(struct Particle *particle_0, struct Particle *particle_1)
{
	if((*particle_0).id >= (*particle_1).id) return;

	double r_01[dim];

	sum_vec(r_01, (*particle_1).pos, (*particle_0).pos, -1);
	double dist = r_abs(r_01);

	if(dist < (*particle_0).radius + (*particle_1).radius)
	{
		r_norm(r_01, r_01);

		double proj0, proj1, vcm;

		proj0 = dot_product((*particle_0).vel, r_01);
		proj1 = dot_product((*particle_1).vel, r_01);

		if(proj0 < proj1) return;

		if(!(*particle_0).fixed && !(*particle_1).fixed)
		{
			vcm = (proj0 * (*particle_0).mass + proj1 * (*particle_1).mass)/
			((*particle_0).mass + (*particle_1).mass);
		}
		else
		{
			vcm = 0;
		}

		sum_vec((*particle_0).vel, (*particle_0).vel, r_01, 2*vcm - 2*proj0);
		sum_vec((*particle_1).vel, (*particle_1).vel, r_01, 2*vcm - 2*proj1);
	}
}
