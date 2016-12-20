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

double zero_vec[dim] = {0};

struct Particle
{
	int id;
	double pos[dim];
	double vel[dim];
	double mass;
	double radius;
	double charge;
	bool fixed;
};

struct point_system
{
	int particleN;
	int forceN;
	int modifierN;
	double time;
	struct Particle Bodies[maxParticles];
	void (*Forces[maxForces])(double*, struct Particle, struct Particle);
	void (*Modifiers[maxModifiers])(struct Particle*, struct Particle*);
};

void
config_system(struct point_system *system)
{
	(*system).time = 0;

	(*system).particleN = 0;
	for (int i = 0; i < maxParticles; ++i)
	{
		if((*system).Bodies[i].mass != 0)
			++(*system).particleN;
	}

	(*system).forceN = 0;
	for (int i = 0; i < maxForces; ++i)
	{
		if((*system).Forces[i] != NULL)
			++(*system).forceN;
	}

	(*system).modifierN = 0;
	for (int i = 0; i < maxModifiers; ++i)
	{
		if((*system).Modifiers[i] != NULL)
			++(*system).modifierN;
	}

	for (int i = 0; i < (*system).particleN; i++)
	{
		(*system).Bodies[i].id = i+1;
	}
}

void
sum_vec(double *c, double *a, double *b, double multiplier_constant)
{
	for(int i = 0; i < dim; i++)
	{
		c[i] = a[i] + multiplier_constant*b[i];
	}
}

double
dot_product(double *a, double *b)
{
	double cnt = 0;
	for(int i = 0; i < dim; i++)
	{
		cnt += a[i]*b[i];
	}
	return cnt;
}

double
r_abs(double *radius)
{
	return sqrt(dot_product(radius, radius));
}

void
r_norm(double *b, double *a)
{
	double norm = r_abs(a);
	for(int i = 0; i < dim; i++)
	{
		b[i] = a[i]/norm;
	}
}

void
calc_system(struct point_system *mysystem, double dt)
{
	for(int i = 0; i < (*mysystem).particleN; i++)
	{
		double force[dim] = {0};

		for(int j = 0; j < (*mysystem).particleN; j++)
		{
			for(int k = 0; k < (*mysystem).forceN; k++)
			{
				(*mysystem).Forces[k](force, (*mysystem).Bodies[i], (*mysystem).Bodies[j]);
			}

			for(int k = 0; k < (*mysystem).modifierN; k++)
			{
				(*mysystem).Modifiers[k](&(*mysystem).Bodies[i], &(*mysystem).Bodies[j]);
			}
		}

		if(!(*mysystem).Bodies[i].fixed)
		{
			sum_vec((*mysystem).Bodies[i].vel, (*mysystem).Bodies[i].vel, force, dt/(*mysystem).Bodies[i].mass);
			sum_vec((*mysystem).Bodies[i].pos, (*mysystem).Bodies[i].pos, (*mysystem).Bodies[i].vel, dt);
		}
	}

	(*mysystem).time += dt;
}

void
print_system(struct point_system system)
{
	printf("%lf", system.time);
	
	for(int i = 0; i < system.particleN; i++)
	{
		for(int j = 0; j < dim; j++)
		{
			printf(" %lf", system.Bodies[i].pos[j]);
		}
	}
	printf("\n");
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
				if(i != j)
				{
					double rad[dim], dist;
					sum_vec(rad, system.Bodies[i].pos, system.Bodies[j].pos, -1);
					dist = r_abs(rad);
					if(dt0*dist/vel < dt)
					{
						dt = dt0*dist/vel;
					}
				}
				
			}
		}
	}
	return fmax(dt, minTimeStep);
}
