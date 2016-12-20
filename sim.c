#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define dim 2
#define maxParticles 10
#define maxForces 3
#define maxModifiers 3

#define maxTimeStep 1e-2
#define minTimeStep 1e-6
#define dt0 1e-4

#include "sim_internal.c"

#define tmax 10

double dynamic_dt(struct point_system);
void grav_force(double *, struct Particle, struct Particle);
void square_wall(struct Particle *, struct Particle *);
void collision(struct Particle *, struct Particle *);
void eletr_force(double *, struct Particle, struct Particle);
void mag_ext_force(double *, struct Particle, struct Particle);

int
main(void)
{
	struct point_system psystem = {
		.Bodies = {
			{
				.pos = {0,0},
				.vel = {1,-3},
				.mass = 50,
				.radius = 0.5,
				.charge = -10,
				.fixed = false
			},
			{
				.pos = {-3,-3},
				.vel = {1,3},
				.mass = 50,
				.radius = 0.5,
				.charge = 10,
				.fixed = false
			},
			{
				.pos = {3,0}, 
				.vel = {-1,-2},
				.mass = 50,
				.radius = 0.5,
				.charge = -10,
				.fixed = false
			},
			{
				.pos = {3,2}, 
				.vel = {-6,0},
				.mass = 50,
				.radius = 0.5,
				.charge = -0,
				.fixed = false
			},
		},
		.Forces = {
			grav_force,
			eletr_force,
			mag_ext_force,
		},
		.Modifiers = {
			square_wall,
			collision,
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

#define grav_c 0.1
void
grav_force(double *force, struct Particle particle_0, struct Particle particle_1)
{
	if(particle_0.id == particle_1.id || particle_0.fixed) return;

	double r_01[dim];
	
	sum_vec(r_01, particle_1.pos, particle_0.pos, -1);

	double dist = r_abs(r_01);
	
	for(int i = 0; i < dim; i++)
	{
		force[i] += + grav_c * particle_0.mass * particle_1.mass * r_01[i] / pow(dist, 3);
	}
}

#define eletr_c 1
void
eletr_force(double *force, struct Particle particle_0, struct Particle particle_1)
{
	if(particle_0.id == particle_1.id || particle_0.fixed) return;

	double r_01[dim];
	
	sum_vec(r_01, particle_1.pos, particle_0.pos, -1);

	double dist = r_abs(r_01);
	
	for(int i = 0; i < dim; i++)
	{
		force[i] += eletr_c * particle_0.charge * particle_1.charge * r_01[i] / pow(dist, 3);
	}
}

#define mag_c 10
void
mag_ext_force(double *force, struct Particle particle_0, struct Particle particle_1)
{
	if(particle_0.id != particle_1.id || particle_0.fixed) return;

	force[0] += mag_c * particle_1.charge * particle_0.vel[1];
	force[1] += - mag_c * particle_1.charge * particle_0.vel[0];
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
collision(struct Particle *particle_0, struct Particle *particle_1)
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
