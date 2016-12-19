#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define dim 2
#define maxParticles 5
#define maxForces 3

#define maxTimeStep 0.01
#define minTimeStep 0.00001
#define dt0 0.0001

#include "sim_internal.c"

#define tmax 10

int main(void)
{
	struct point_system psystem = {
		.Bodies = {
			{.pos = {5.3,0}, .vel = {0,1}, .mass = 2, .fixed = false},
			{.pos = {4.7,0}, .vel = {0,-1}, .mass = 2, .fixed = false},
			{.pos = {0,0}, .vel = {0,-2}, .mass = 100, .fixed = false},
			{.pos = {10,0}, .vel = {0,2}, .mass = 100, .fixed = false},
		},
		.Forces = {
			grav_force,
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
