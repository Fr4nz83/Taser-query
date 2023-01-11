#pragma once


#include "UndirectedRN.hpp"
#include "lemon/list_graph.h"
#include "lemon/smart_graph.h"


template<typename GR>
class EstimateCalculator
{
protected:

	// *** PROTECTED TYPES *** //

	typedef GR Graph;
	typedef typename Graph::Node Node;
	typedef typename Graph::template NodeMap<typename std::pair<double,double>> NodeMap;
	typedef typename UndirectedRN::Spatial2DPoint Spatial2DPoint;


	// *** PROTECTED FIELDS *** //

	static const GeographicLib::Geodesic geod;
	NodeMap& mapVertices;



	// *** PROTECTED METHODS *** //

	inline double computeDistancePoints(const Spatial2DPoint& p1, const Spatial2DPoint& p2) const
	{
		// static const GeographicLib::Geodesic& geod(GeographicLib::Geodesic::WGS84());
		double distance_metres = 0;
		EstimateCalculator::geod.Inverse(p1.first, p1.second,
					 	 	 	 	 	 p2.first, p2.second, distance_metres);
		return(distance_metres);
	}



public:

	// *** PUBLIC CTORS/DTOR *** //

	EstimateCalculator(NodeMap& mapVertices) : mapVertices(mapVertices) {};



	// *** PUBLIC METHODS *** //

	inline double operator()(const Node& u, const Node& v) const
	{
		return this->computeDistancePoints(this->mapVertices[u], this->mapVertices[v]);
	}
};


// *** INITIALIZATION STATIC FIELDS CLASS *** //

template<typename GR>
const GeographicLib::Geodesic EstimateCalculator<GR>::geod = GeographicLib::Geodesic::WGS84();
