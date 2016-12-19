struct Particle
{
	int id;
	double pos[dim];
	double vel[dim];
	double mass;
	bool fixed;
};

struct point_system
{
	int particleN;
	int forceN;
	double time;
	struct Particle Bodies[maxParticles];
	void (*Forces[maxForces])(double*, struct Particle, struct Particle);
};

void initSystem(struct point_system *system)
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

	for (int i = 0; i < (*system).particleN; i++)
	{
		(*system).Bodies[i].id = i;
	}
}

void sum_vec(double *c, double *a, double *b, double multiplier_constant)
{
	for(int i = 0; i < dim; i++)
	{
		c[i] = a[i] + multiplier_constant*b[i];
	}
}

double r_abs(double *radius)
{
	double norm = 0;
	for(int i = 0; i < dim; i++)
	{
		norm += radius[i]*radius[i];
	}
	return sqrt(norm);
}

#define grav_c 1
void grav_force(double *force, struct Particle particle_0, struct Particle particle_1)
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

void calc_system(struct point_system *mysystem, double dt)
{
	for(int i = 0; i < (*mysystem).particleN; i++)
	{
		if((*mysystem).Bodies[i].fixed)
		{
			i++;
			continue;
		}
		
		double force[dim] = {0};
		
		for(int j = 0; j < (*mysystem).particleN; j++)
		{
			for(int k = 0; k < (*mysystem).forceN; k++)
			{
				(*mysystem).Forces[k](force, (*mysystem).Bodies[i], (*mysystem).Bodies[j]);
			}
		}
		
		sum_vec((*mysystem).Bodies[i].vel, (*mysystem).Bodies[i].vel, force, dt/(*mysystem).Bodies[i].mass);
		sum_vec((*mysystem).Bodies[i].pos, (*mysystem).Bodies[i].pos, (*mysystem).Bodies[i].vel, dt);
	}

	(*mysystem).time += dt;
}

void print_system(struct point_system system)
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

// Idea:
// dt = dt0 * min(dist/vel) where dt0 << 1, so that vel*dt << dist
#define maxTimeStep 0.01
#define minTimeStep 0.00001
#define dt0 0.0001
double dynamic_dt(struct point_system system)
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
