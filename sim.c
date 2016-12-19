#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define dim 2
#define maxParticles 5
#define maxForces 3
#define maxModifiers 3

#define maxTimeStep 1e-2
#define minTimeStep 1e-6
#define dt0 1e-4

#include "sim_internal.c"

#define tmax 30

double dynamic_dt(struct point_system);
void grav_force(double *, struct Particle, struct Particle);
void square_wall(struct Particle *, struct Particle *);
void collision(struct Particle *, struct Particle *);

int
main(void)
{
	struct point_system psystem = {
		.Bodies = {
			{.pos = {0,0}, .vel = {0,1}, .mass = 80, .radius = .5, .fixed = false},
			{.pos = {3,0}, .vel = {-1,-2}, .mass = 40, .radius = .5, .fixed = false},
			{.pos = {6,0}, .vel = {3,0}, .mass = 100, .radius = .3, .fixed = false},
			{.pos = {9,0}, .vel = {-2,3}, .mass = 200, .radius = .2, .fixed = false},
		},
		.Forces = {
			grav_force,
		},
		.Modifiers = {
			square_wall,
			collision,
		},
	};

	initSystem(&psystem);

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
	if(particle_0.id == particle_1.id) return;

	double r_01[dim];
	
	sum_vec(r_01, particle_1.pos, particle_0.pos, -1);

	double dist = r_abs(r_01);
	
	for(int i = 0; i < dim; i++)
	{
		force[i] = force[i] + grav_c * particle_0.mass * particle_1.mass * r_01[i] / pow(dist, 3);
	}
}

#define topw 5
#define bottomw -5
#define rightw 12
#define leftw -2
void 
square_wall(struct Particle *particle_0, struct Particle *particle_1)
{
	if((*particle_0).id != (*particle_1).id) return;

	if((*particle_0).pos[1] > topw && (*particle_0).vel[1] > 0)
		(*particle_0).vel[1] *= -1;

	if((*particle_0).pos[0] > rightw && (*particle_0).vel[0] > 0)
		(*particle_0).vel[0] *= -1;

	if((*particle_0).pos[1] < bottomw && (*particle_0).vel[1] < 0)
		(*particle_0).vel[1] *= -1;

	if((*particle_0).pos[0] < leftw && (*particle_0).vel[0] < 0)
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

		sum_vec((*particle_0).vel, (*particle_0).vel, r_01, -proj0);
		sum_vec((*particle_1).vel, (*particle_1).vel, r_01, -proj1);

		vcm = (proj0 * (*particle_0).mass + proj1 * (*particle_1).mass)/
		((*particle_0).mass + (*particle_1).mass);
		proj0 = 2*vcm - proj0;
		proj1 = 2*vcm - proj1;

		sum_vec((*particle_0).vel, (*particle_0).vel, r_01, proj0);
		sum_vec((*particle_1).vel, (*particle_1).vel, r_01, proj1);
	}
}

// dt = dt0 * min(dist/vel) where dt0 << 1, so that vel*dt << dist
double
dynamic_dt(struct point_system system)
{
	double dt = maxTimeStep;

	for (int i = 0; i < system.particleN; i++)
	{
		if(system.Bodies[i].fixed == false)
		{
			double vel = r_abs(system.Bodies[i].vel);
			for (int j = 0; j < system.particleN; j++)
			{
				if(i == j) break;
				double rad[dim], dist;
				sum_vec(rad, system.Bodies[i].pos, system.Bodies[j].pos, -1);
				dist = r_abs(rad);
				if(dt0*dist/vel < dt)
					dt = dt0*dist/vel;
			}
		}
	}
	return fmax(dt, minTimeStep);
}
