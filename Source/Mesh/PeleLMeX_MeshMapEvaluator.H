#ifndef PELELMEX_MESHMAP_EVALUATOR_H
#define PELELMEX_MESHMAP_EVALUATOR_H

// The mesh-mapping descriptor lives in PelePhysics
// (Source/Utility/MeshMap/MeshMapEvaluator.H) so that the turbulent-inflow
// reader and the inflow-file generators share the solver's definition of
// the mapping.  PeleLMeX code keeps referring to it by its unqualified name.
#include <MeshMapEvaluator.H>

using MeshMapEvaluator = pele::physics::MeshMapEvaluator;

#endif // PELELMEX_MESHMAP_EVALUATOR_H
