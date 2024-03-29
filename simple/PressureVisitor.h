#ifndef SIMPLE_PRESSUREVISITOR_H_
#define SIMPLE_PRESSUREVISITOR_H_

#include "simple.decl.h"
#include "common.h"
#include <cmath>

struct PressureVisitor {
// in leaf check for not same particle plz
private:
  const Real radius = .01;

public:
  bool node(SourceNode<CentroidData> source, TargetNode<CentroidData> target) {
    Real rsq = (source.data->getCentroid() - target.data->getCentroid()).lengthSquared();
    return (rsq < radius * radius);
    // this just looks at the centroids when it should look at the whole boxes
    // need to use intersect(), maybe use Sphere<> ?
  }

  void leaf(SourceNode<CentroidData> source, TargetNode<CentroidData> target) {
    for (int i = 0; i < target.n_particles; i++) {
      for (int j = 0; j < source.n_particles; j++) {
        if ((target.particles[i].position - source.particles[j].position).lengthSquared() < radius * radius) {
          //import kernel math here
        }
      }
    }
  }
};

#endif // SIMPLE_PRESSUREVISITOR_H_
