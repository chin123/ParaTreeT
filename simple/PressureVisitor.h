#ifndef SIMPLE_PRESSUREVISITOR_H_
#define SIMPLE_PRESSUREVISITOR_H_

#include "simple.decl.h"
#include "common.h"
#include <cmath>
#include <iostream>

struct PressureVisitor {
// in leaf check for not same particle plz
private:
  const Real radius = 100;
  SPH::SplineKernel sp;
  const Real restDensity = 1000.0; // Should change based on what substance is being simulated
  const Real gasConstant = 2000.0; // Should change based on temperature

  bool intersect(SourceNode<CentroidData>& source, TargetNode<CentroidData>& target) {
    Real rsq = radius * radius;
    Vector3D<Real> pos = source.data->getCentroid();
    const OrientedBox<Real> box = target.data->box;
    double dsq = 0.0;
    double delta;

    if((delta = box.lesser_corner.x - pos.x) > 0)
	dsq += delta * delta;
    else if((delta = pos.x - box.greater_corner.x) > 0)
	dsq += delta * delta;
    if(rsq < dsq)
	return false;
    if((delta = box.lesser_corner.y - pos.y) > 0)
	dsq += delta * delta;
    else if((delta = pos.y - box.greater_corner.y) > 0)
	dsq += delta * delta;
    if(rsq < dsq)
	return false;
    if((delta = box.lesser_corner.z - pos.z) > 0)
	dsq += delta * delta;
    else if((delta = pos.z - box.greater_corner.z) > 0)
	dsq += delta * delta;
    //std::cout << "dsq: " << dsq << " rsq: " << rsq << std::endl;
    return (dsq <= rsq);
  }

public:
  bool node(SourceNode<CentroidData> source, TargetNode<CentroidData> target) {
    //Real rsq = (source.data->getCentroid() - target.data->getCentroid()).lengthSquared();
    // return (rsq < radius * radius);
    // this just looks at the centroids when it should look at the whole boxes
    // need to use intersect(), maybe use Sphere<> ?

    //std::cout << "intersect? " << intersect(source, target) << std::endl;
    //return intersect(source, target);
    return true;
  }

  void leaf(SourceNode<CentroidData> source, TargetNode<CentroidData> target) {
    std::cout << "Reached leaf" << std::endl;
    for (int i = 0; i < target.n_particles; i++) {
      Vector3D<Real> force(0, 0, 0);
      Real pi = gasConstant * (target.particles[i].density - restDensity);
      for (int j = 0; j < source.n_particles; j++) {
        if (target.particles[i].key == source.particles[j].key) continue;
        std::cout << "example radius " << (target.particles[i].position - source.particles[j].position).lengthSquared() << std::endl;
        if ((target.particles[i].position - source.particles[j].position).lengthSquared() < radius * radius) {
          //import kernel math here
          Real pj = gasConstant * (source.particles[j].density - restDensity);
          Real gradient = sp.evaluateGradient(sqrt((target.particles[i].position - source.particles[j].position).lengthSquared()), radius);
          force += source.particles[j].mass * ((pi + pj)/(2*source.particles[j].density)) * gradient;

        }
      }
      force *= -1;
      target.applyForce(i, force);
    }
  }
};

#endif // SIMPLE_PRESSUREVISITOR_H_
