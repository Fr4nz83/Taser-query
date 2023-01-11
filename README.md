In order to correctly compile and execute the application, the following shared libraries need to be present:
- Boost
- Intel TBB
- GeographicLib

Example of a single run: 

./bin/CASRQ -v datasets/Oslo/roadnetwork/RoadVerticesOSLO.txt -e datasets/Oslo/roadnetwork/RoadEdgesOSLO.txt -p datasets/Oslo/poi/uniform/AtmBanksOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/CoffeeShopsOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/GasStationsOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/MovieTheatersOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/PharmaciesOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/PubsBarsOSLO.txt.CORRECTED.txt.UNIFORM.txt datasets/Oslo/poi/uniform/RestaurantsOSLO.txt.CORRECTED.txt.UNIFORM.txt -a 0 -m 2 -r 10
