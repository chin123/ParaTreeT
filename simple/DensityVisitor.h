#ifndef SIMPLE_DENSITYVISITOR_H_
#define SIMPLE_DENSITYVISITOR_H_

#include "simple.decl.h"
#include "common.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <queue>
//#include "ParticleComp.h"

struct DensityVisitor {
// in leaf check for not same particle plz
private:
  const int k = 32;
  const Real radius = 100;
  SPH::SplineKernel sp;
private:

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
    return (dsq <= rsq);
  }

public:
  bool node(SourceNode<CentroidData> source, TargetNode<CentroidData> target) {
    return intersect(source, target);
  }

  void leaf(SourceNode<CentroidData> source, TargetNode<CentroidData> target) {
    for (int i = 0; i < target.n_particles; i++) {
      double density = 0;
      for (int j = 0; j < source.n_particles; j++) {
        Vector3D<Real> diff = target.particles[i].position - source.particles[j].position;
        if (diff.lengthSquared() <= radius * radius) {
          Vector3D<Real> diff = source.particles[j].position - target.particles[i].position;
          density += source.particles[j].mass * sp.evaluate(sqrt(diff.lengthSquared()), radius);
        }
      }
      target.particles[i].density += density;
    }
  }
};

#endif // SIMPLE_DENSITYVISITOR_H_
