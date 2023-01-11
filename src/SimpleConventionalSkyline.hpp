#include "SimpleSkyline.hpp"


class SimpleConventionalSkyline : public SimpleSkyline
{
public:

	// *** PUBLIC CTORS/DTORS *** //

	SimpleConventionalSkyline(const ElementSkyline& ub1, const std::vector<uint32_t>& ub1Seq);
	~SimpleConventionalSkyline() {};



	// *** PUBLIC METHODS *** //

	virtual bool insertElement(const ElementSkyline& el, const std::vector<uint32_t>& elSeq = std::vector<uint32_t>(0));
};
