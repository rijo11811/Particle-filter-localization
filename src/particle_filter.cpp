/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <cassert>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  num_particles = 100;  // TODO: Set the number of particles
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  for (int i = 0; i < num_particles; ++i) {
    Particle p;
    p.id = i;
    p.x =dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight =1.0;
    particles.push_back(p);
    particles[i].weight = 1.0;
  } 
  is_initialized = true;
    //std::cout<<particles[0].x<<" "<<particles[0].y<<" "<<particles[0].theta<<" "<<particles[0].weight<<"\n";

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine gen;
    
	for (int i = 0; i < num_particles; ++i){
      
        if (fabs(yaw_rate) < 0.0001){
            particles[i].x += velocity*delta_t * cos(particles[i].theta);
            particles[i].y += velocity*delta_t * sin(particles[i].theta);
            
        }
        else{
            
            particles[i].x += velocity/yaw_rate * (sin(particles[i].theta +  delta_t * yaw_rate) - sin(particles[i].theta));
            particles[i].y += velocity/yaw_rate * (-cos(particles[i].theta +  delta_t * yaw_rate) + cos(particles[i].theta));
            particles[i].theta += yaw_rate * delta_t;
          
          
        }
        
      
        // Add random Gaussian noise
      std::normal_distribution<double> dist_x(particles[i].x , std_pos[0]);
	  std::normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
	  std::normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
      particles[i].x = dist_x(gen);
      particles[i].y = dist_y(gen);
      particles[i].theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for(unsigned int i = 0; i < observations.size(); ++i)
  {
    double min_dist = std::numeric_limits<double>::max();
    int position = -1;

    for(unsigned int j = 0; j < predicted.size(); ++j)
    { double distance = 0.0;
      distance= dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
      if(distance <  min_dist)
      {
        min_dist = distance;
        position = j;
      }
        
    }
    observations[i].id=predicted[position].id;
    //std::cout<<"asso"<<observations[i].x<<" "<<observations[i].y<<" "<<predicted[position].x<<" "<<predicted[position].y<<"\n";
    
        
    }
  

}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   
   vector<int> associations;
   vector<double> sense_x;
   vector<double> sense_y;
   std::vector<LandmarkObs> landmark_predictions;
   double gauss_norm;
   gauss_norm = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
   for( int i = 0; i <num_particles; ++i)
   { 
     double x_part, y_part;
     x_part = particles[i].x;
     y_part = particles[i].y;
     landmark_predictions.clear();
     for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j)
     {
       if(dist(x_part,y_part,map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f) <=sensor_range)
       {  
         struct LandmarkObs l;
         l.id = map_landmarks.landmark_list[j].id_i;
         l.x=static_cast<double>(map_landmarks.landmark_list[j].x_f);
         l.y=static_cast<double>(map_landmarks.landmark_list[j].y_f);
         landmark_predictions.push_back(l);
       }
    }
     assert(!landmark_predictions.empty());
    
    std::vector<LandmarkObs>temp_obs;
    temp_obs.clear();
    for (unsigned int j = 0; j <observations.size(); ++j)
    { LandmarkObs l;
      l.x=x_part + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
      l.y=y_part + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
      
      temp_obs.push_back(l);
        
    }
    dataAssociation(landmark_predictions,temp_obs);   
    //SetAssociations(particles[i],temp_obs[],const vector<double>& sense_x,const vector<double>& sense_y)
     /*
     if(i==10)
     {std::cout<<" lm "<<landmark_predictions.size()<<"obs s"<<temp_obs.size()<<std::endl;}
     */
     particles[i].weight = 1.0;
     for(unsigned int j = 0; j < temp_obs.size(); j++)
     { double w=0.0;
       double x_map = 0.0;
  	   x_map = temp_obs[j].x;
 	   double y_map = 0.0;
  	   y_map = temp_obs[j].y;    
  	   double exponent=0.0;
       double mu_x=0.0,mu_y=0.0;
      /*
       for(unsigned int k = 0; k < landmark_predictions.size(); ++k)
       {
         if(temp_obs[j].id == landmark_predictions[k].id)
         {
           mu_x = landmark_predictions[k].x;
           mu_y = landmark_predictions[k].y;
           //std::cout<<"mu";
         }
         
       }
      */
      mu_x= static_cast<double>(map_landmarks.landmark_list[temp_obs[j].id-1].x_f);
      mu_y=static_cast<double>(map_landmarks.landmark_list[temp_obs[j].id-1].y_f);
       double std_x_2 = pow(std_landmark[0], 2.0);
       double std_y_2 = pow(std_landmark[1], 2.0);
       double dx,dy;
       dx = x_map - mu_x;
       dy = y_map - mu_y;
  	   exponent = (pow(dx, 2) / (2 * std_x_2)) + (pow(dy, 2) / (2 * std_y_2));     		
      //std::cout<<gauss_norm<<"  "<<exponent<<" "<<exp(-exponent)<<" "<<x_map<<" "<<mu_x<<" "<<y_map<<" "<<mu_y<<"\n";
      /*
      if(i==10)
      {
        std::cout<<"x_map - mu_x"<<x_map - mu_x<<"  "<<"y_map - mu_y"<<y_map - mu_y<<" "<<"exponent"<<exponent<<std::endl;
      }
      */
       w = gauss_norm * exp(-1*exponent);
      
      
      particles[i].weight *=  w;
      
       //std::cout<<map_landmarks.landmark_list[temp_obs[j].id-1].x_f<<"\n";
       //std::cout<<mu_x<<" "<<mu_y<<" "<<temp_obs[3].id<<" "<<temp_obs.size()<<"\n";
          
      	associations.push_back(temp_obs[j].id);
		sense_x.push_back(temp_obs[j].x);
		sense_y.push_back(temp_obs[j].y);
     }
     weights.push_back(particles[i].weight);
     SetAssociations(particles[i], associations, sense_x, sense_y);
     associations.clear();
     sense_x.clear();
     sense_y.clear();
     
   }
  

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
    std::random_device rd;
    std::mt19937 gen(rd()); 	
  	std::discrete_distribution<> d(weights.begin(),weights.end());
    std::vector<Particle> particles_temp;
    particles_temp.clear();
    particles_temp.resize(num_particles);
    for( int i = 0; i <num_particles; ++i)
    {
      //particles_temp.push_back(particles[d(gen)]);
      particles_temp[i] = particles[d(gen)];
    }
    
    particles = particles_temp;
  	//std::cout<<"w size"<<weights.size()<<"W[100]"<<weights[100]<<std::endl;
  	weights.clear();

 
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}