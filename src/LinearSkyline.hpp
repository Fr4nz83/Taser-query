/**
 * @brief Class modeling a linear skyline.
 * @author Francesco Lettich
 */

#pragma once


// *** INCLUDES *** //

#include <set>
#include <iostream>
#include <unordered_map>
#include <vector>


class LinearSkyline
{
public:

	// *** PUBLIC TYPEDEFS *** //

	typedef typename std::pair<double,double> ElementSkyline;
	template<typename K>
	struct LSCmp {inline bool operator()(const K& lhs, const K& rhs) {return lhs.first < rhs.first;}};
	typedef std::set<ElementSkyline, LSCmp<ElementSkyline>> Skyline;



protected:

	// *** PROTECTED FIELDS *** //

	Skyline LS;    // Sorted list containing the linear skyline computed so far.
	std::unordered_map<double, std::vector<uint32_t>> tableSequences; // Lookup table containing the sequences
																	  // associated with the elements in the skyline.



public:

	// *** PUBLIC CTORS/DTORS *** //

	LinearSkyline(){};
	LinearSkyline(const ElementSkyline& ub1, const std::vector<uint32_t>& ub1Seq,
				  const ElementSkyline& ub2, const std::vector<uint32_t>& ub2Seq);
	~LinearSkyline() {};



	// *** PUBLIC METHODS *** //

	Skyline::const_iterator findRightNeighbor(const ElementSkyline& el) {return(this->LS.upper_bound(el));};
	Skyline::const_iterator cbegin() {return(this->LS.cbegin());}
	Skyline::const_iterator cend() {return(this->LS.cend());}
	inline const ElementSkyline& getFirstUpperBound() {return(*(this->LS.begin()));}
	inline const ElementSkyline& getSecondUpperBound() {return(*(this->LS.rbegin()));}
	bool insertElement(const ElementSkyline& el, const std::vector<uint32_t>& elSeq = std::vector<uint32_t>(0));
	const Skyline& getLS() const {return(this->LS);};
	const std::vector<uint32_t>& getSequence(const double& POICost) {return this->tableSequences[POICost];};
	void mergeSkyline(LinearSkyline& LSOther);

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
inline bool LinearSkyline::checkDomination(const ElementSkyline& l1, const ElementSkyline& l2)
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
inline bool LinearSkyline::checkNormalDomination(const ElementSkyline& l1, const ElementSkyline& l2, const ElementSkyline& p)
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
