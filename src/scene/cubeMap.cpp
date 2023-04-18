#include "cubeMap.h"
#include "ray.h"
#include "../ui/TraceUI.h"
#include "../scene/material.h"
extern TraceUI* traceUI;

glm::dvec3 CubeMap::getColor(ray r) const
{
	// YOUR CODE HERE
	// FIXME: Implement Cube Map here
	glm::dvec3 direction = r.getDirection();
	if (glm::abs(direction.x) > glm::abs(direction.y) && glm::abs(direction.x) > glm::abs(direction.z)) {
		double y = 0.5 * (direction.y/glm::abs(direction.x) + 1);
		double z = 0.5 * (direction.z/glm::abs(direction.x) + 1);
		if (direction.x > 0) {
			glm::dvec2 coords(z, y);
			return tMap[0]->getMappedValue(coords);
		}
		else {
			glm::dvec2 coords(1 - z, y);
			return tMap[1]->getMappedValue(coords);
		}
	}
	if (glm::abs(direction.y) > glm::abs(direction.x) && glm::abs(direction.y) > glm::abs(direction.z)) {
		double x = 0.5 * (direction.x/glm::abs(direction.y) + 1);
		double z = 0.5 * (direction.z/glm::abs(direction.y) + 1);
		if (direction.y > 0) {
			glm::dvec2 coords(x, z);
			return tMap[2]->getMappedValue(coords);
		}
		else {
			glm::dvec2 coords(x, 1 - z);
			return tMap[3]->getMappedValue(coords);
		}
	}
	if (glm::abs(direction.z) > glm::abs(direction.x) && glm::abs(direction.z) > glm::abs(direction.y)) {
		double x = 0.5 * (direction.x/glm::abs(direction.z) + 1);
		double y = 0.5 * (direction.y/glm::abs(direction.z) + 1);
		if (direction.z < 0) {
			glm::dvec2 coords(x, y);
			return tMap[4]->getMappedValue(coords);
		}
		else {
			glm::dvec2 coords(1 - x, y);
			return tMap[5]->getMappedValue(coords);
		}
	}
	return glm::dvec3();
}

CubeMap::CubeMap()
{
}

CubeMap::~CubeMap()
{
}

void CubeMap::setNthMap(int n, TextureMap* m)
{
	if (m != tMap[n].get())
		tMap[n].reset(m);
}
