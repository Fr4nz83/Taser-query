/**
 * @brief Class modeling a linear skyline.
 * @author Francesco Lettich
 */

#pragma once


// *** INCLUDES *** //

#include <vector>
#include <iostream>


class SimpleSkyline
{
public:

	// *** PUBLIC TYPEDEFS *** //

	typedef typename std::pair<double,double> ElementSkyline;
	typedef std::vector<ElementSkyline> Skyline;
	typedef std::vector<std::vector<uint32_t>> SkylineSequences;



protected:

	// *** PROTECTED FIELDS *** //

	Skyline LS;    			 // Linear skyline.
	SkylineSequences LSS;    // Vector containing the sequences associated with the elements of the lin. skylines.



	// *** PROTECTED CTORS *** //

	SimpleSkyline(){}; // Default ctor, not callable from outside!



public:

	// *** PUBLIC CTORS/DTORS *** //

	SimpleSkyline(const ElementSkyline& ub1);
	SimpleSkyline(const ElementSkyline& ub1, const std::vector<uint32_t>& ub1Seq);
	~SimpleSkyline() {};



	// *** PUBLIC METHODS *** //

	Skyline::const_iterator cbegin() {return(this->LS.cbegin());}
	Skyline::const_iterator cend() {return(this->LS.cend());}
	virtual bool insertElement(const ElementSkyline& el, const std::vector<uint32_t>& elSeq = std::vector<uint32_t>(0));
	Skyline getLS() {return(this->LS);};
	SkylineSequences getLSS() {return(this->LSS);};
	SimpleSkyline mergeWithSkyline(const SimpleSkyline& SS);

	static bool checkDomination(const ElementSkyline& l1, const ElementSkyline& l2);
	static bool checkNormalDomination(const ElementSkyline& l1, const ElementSkyline& l2, const ElementSkyline& p);

	// Internal test methods.
	static void testLinSky();
};



//*** INLINE PUBLIC METHODS DECLARATIONS ***//

/**
 * @brief This static method takes in input two bidimensional points and checks whether l1 dominates l2.
 *
 * @return "true" if l1 dominates l2, "false" otherwise.
 */
inline bool SimpleSkyline::checkDomination(const ElementSkyline& l1, const ElementSkyline& l2)
{
	return(((l1.first <= l2.first) && (l1.second < l2.second)) ||
		   ((l1.first < l2.first) && (l1.second <= l2.second)));
}

/**
 * @brief This static method takes in input two points and checks whether l1 and l2 linearly dominate p.
 *
 * @note with respect to the X axis, l1 has to be the "left" neighbor of p, while l2 has to be the "right" neighbor.
 * 		 As such, the following conditions must hold: l2.x >= l1.x AND l1.y >= l2.y.
 */
inline bool SimpleSkyline::checkNormalDomination(const ElementSkyline& l1, const ElementSkyline& l2, const ElementSkyline& p)
{
	// Component-wise minimum of "l1" and "l2".
	const ElementSkyline u = std::make_pair(std::min(l1.first, l2.first), std::min(l1.second, l2.second));

	// Compute the normal of the segment connecting l1 and l2.
	const ElementSkyline normal = std::make_pair(-(l2.second - l1.second), (l2.first - l1.first));

//	std::cout << "Component-wise minimum: (" << u.first << "," << u.second << ")\n";
//	std::cout << "LinComb el: " << (normal.first * p.first + normal.second * p.second) << "\n";
//	std::cout << "LinComb sky: " << (normal.first * l1.first + normal.second * l1.second) << "\n";

	// Check linear dominance of l1 and l2 over p.
	const bool test = ((u.first < p.first && u.second <= p.second) || (u.first <= p.first && u.second < p.second)) &&
				((normal.first * p.first + normal.second * p.second) >
				 (normal.first * l1.first + normal.second * l1.second));

	return(test);
}
