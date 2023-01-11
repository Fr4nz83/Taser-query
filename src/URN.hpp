#pragma once


// *** INCLUDES *** //
#include "UndirectedRN.hpp"



class URN : public UndirectedRN
{
protected:

	// *** PROTECTED TYPES *** //

	typedef struct
	{
		uint32_t idEdge;
		double distStart;
		double cost;
		uint32_t idOrgVertex;
	} POIFile;



	// *** PROTECTED FIELDS *** //

	// Backup original graph.
	Graph gBackup;
	std::vector<Graph::Node> vecNodesBackup;
	Graph::NodeMap<Spatial2DPoint> mapVerticesBackup;
	std::vector<Graph::Edge> vecEdgesBackup;
	Graph::EdgeMap<double> mapEdgesBackup;

	// Backup initial TD hash map.
	std::unordered_map<uint32_t, std::vector<std::pair<double, uint32_t>>> tableTDsBackup;

	// Vectors used to extract random POIs from the given COIs.
	std::vector<std::vector<POIFile>> vecCOIsBackup;
	std::vector<std::vector<POIFile>> vecCOIs;



	// *** PROTECTED METHODS *** //

	void copyGraph(Graph& gOld, std::vector<Graph::Node>& vecNodesOld, Graph::NodeMap<Spatial2DPoint>& mapVerticesOld, std::vector<Graph::Edge>& vecEdgesOld, Graph::EdgeMap<double>& mapEdgesOld,
				   Graph& gNew, std::vector<Graph::Node>& vecNodesNew, Graph::NodeMap<Spatial2DPoint>& mapVerticesNew, std::vector<Graph::Edge>& vecEdgesNew, Graph::EdgeMap<double>& mapEdgesNew);

	void parseFilePOI(const std::vector<const char*>& vecNameFilesPOIs);

	void applyPOIs(uint32_t numPOIsCOI);
	void updateGraph();



public:

	// *** PUBLIC CTORS/DTOR *** //

	URN(const char* nameFileVertices,
	 	const char* nameFileEdges,
		std::vector<const char*>& vecNameFilesPOIs);
	~URN() {};



	// *** PUBLIC METHODS *** //

	void loadPrecomputedTDPOIs(const std::string& nameFile);
	void regenerateGraph(const uint32_t& numPOIsCOI);
};
