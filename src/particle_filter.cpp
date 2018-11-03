/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <limits>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	// create normal distribution for sampling
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// set num of particles
	num_particles = 10;

	// create all particles
	for (int i=0; i < num_particles; ++i) {
		// create single particle
		Particle particle;	
		particle.id = i+1;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);

		// init weight vector
		weights.push_back(1);		
	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	// random generator for gaussian noise
	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	// check for 0 yaw rate
	if (yaw_rate != 0) {
		for (int i=0; i < num_particles; ++i) {
			// intermediate variables
			double theta = particles[i].theta;
			double v_yaw = velocity / yaw_rate;
			
			// calculate new particle state
			double nx = particles[i].x + v_yaw*(sin(theta + yaw_rate*delta_t) - sin(theta));
			double ny = particles[i].y + v_yaw*(cos(theta) - cos(theta + yaw_rate*delta_t));
			double ntheta = particles[i].theta + yaw_rate*delta_t;

			// add noise
			nx += dist_x(gen);
			ny += dist_y(gen);
			ntheta += dist_theta(gen);	

			// set particle position
			particles[i].x = nx; 
			particles[i].y = ny;
			particles[i].theta = ntheta; 
		}
	}
	else {
		for (int i=0; i < num_particles; ++i) {
			// calculate new particle state
			double ntheta = particles[i].theta;
			double nx = particles[i].x + velocity*delta_t*cos(ntheta);
			double ny = particles[i].y + velocity*delta_t*sin(ntheta);
			
			// sample new state from gaussian
			normal_distribution<double> dist_nx(nx, std_pos[0]);
			normal_distribution<double> dist_ny(ny, std_pos[1]);
			normal_distribution<double> dist_ntheta(ntheta, std_pos[2]);

			particles[i].x = dist_nx(gen);
			particles[i].y = dist_ny(gen);
			particles[i].theta = dist_ntheta(gen);
		}	
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	// gaussian noise for observation coordinates
	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_landmark[0]);
	normal_distribution<double> dist_y(0, std_landmark[1]);

	// loop through each particle to update
	for (int i=0; i < num_particles; ++i) {
		for (unsigned int j=0; j < observations.size(); ++j) {
			/////// TRANSFORM //////
			LandmarkObs map_coords;
			// intermediate variables
			double x_particle = particles[i].x;
			double y_particle = particles[i].y;
			double x_obs = observations[j].x + dist_x(gen);
			double y_obs = observations[j].y + dist_y(gen);
			double theta = particles[i].theta;

			// calculate map coords
			map_coords.x = double(x_particle + (cos(theta)*x_obs) - (sin(theta)*y_obs));
			map_coords.y = double(y_particle + (sin(theta)*x_obs) + (cos(theta)*y_obs));

			/////// ASSOCIATE /////////
			// find closest landmark for this transformed observation 
			double shortest_distance = std::numeric_limits<double>::infinity(); 
			for (unsigned int k=0; k < map_landmarks.landmark_list.size(); ++k) {
				// intermediate variables
				double land_x = map_landmarks.landmark_list[k].x_f;
				double land_y = map_landmarks.landmark_list[k].y_f;
				double distance = dist(map_coords.x, map_coords.y, land_x, land_y);
				
				if (distance < shortest_distance) {
					shortest_distance = distance;
					map_coords.id = map_landmarks.landmark_list[k].id_i;
				}
			}
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
