#include <cmath>
#include <iostream>

#include "light.h"
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>


using namespace std;

double DirectionalLight::distanceAttenuation(const glm::dvec3& P) const
{
	// distance to light is infinite, so f(di) goes to 0.  Return 1.
	printf("test\n");
	return 1.0;
}


glm::dvec3 DirectionalLight::shadowAttenuation(const ray& r, const glm::dvec3& p) const
{
	printf("directionallight shadow atten\n");
	// YOUR CODE HERE:
	glm::dvec3 attenuation = r.getAtten();
	glm::dvec3 rayDirection = r.getDirection();
	glm::dvec3 directionOfLightFromRayIntersect = this->getDirection(p);
	ray shadow(p + RAY_EPSILON * (-rayDirection), directionOfLightFromRayIntersect, glm::dvec3(1,1,1), ray::SHADOW);
	isect shadowIntersect;
	bool intersecting = scene->intersect(shadow, shadowIntersect);
	glm::dvec3 shadowIntersectPosition = shadow.at(shadowIntersect);
	glm::dvec3 directionOfLightFromShadowIntersect = getDirection(shadowIntersectPosition);
	glm::dvec3 previousRayDirection = rayDirection;
	while (intersecting) {
		Material m = shadowIntersect.getMaterial();
		if (m.Trans()) {
			bool leaving = glm::dot(directionOfLightFromShadowIntersect, shadowIntersect.getN()) > 0;
			if (leaving) {
				double d = glm::distance(shadow.getPosition(), shadowIntersectPosition);
				attenuation *= glm::pow(m.kt(shadowIntersect), glm::dvec3(d));
			}
		}
		else {
			attenuation = glm::dvec3(0.0, 0.0, 0.0);
			break;
		}
		previousRayDirection = directionOfLightFromShadowIntersect;
		directionOfLightFromShadowIntersect = getDirection(shadowIntersectPosition);
		shadow = ray(shadowIntersectPosition + RAY_EPSILON * previousRayDirection, directionOfLightFromShadowIntersect, glm::dvec3(1,1,1), ray::SHADOW);
		intersecting = scene->intersect(shadow, shadowIntersect);
		shadowIntersectPosition = shadow.at(shadowIntersect);
	}
	printf("directionallight shadow atten finish\n");
	return attenuation;
}

glm::dvec3 DirectionalLight::getColor() const
{
	printf("color\n");
	return color;
}

glm::dvec3 DirectionalLight::getDirection(const glm::dvec3& P) const
{
	return -orientation;
}

double PointLight::distanceAttenuation(const glm::dvec3& P) const
{

	// YOUR CODE HERE

	// You'll need to modify this method to attenuate the intensity 
	// of the light based on the distance between the source and the 
	// point P.  For now, we assume no attenuation and just return 1.0
	// f(d) = min( 1, 1/( a + b d + c d^2 ) )
	// float constantTerm;		// a
	// float linearTerm;		// b
	// float quadraticTerm;		// c
	printf("pointlight dist atten\n");
	double d = sqrt(pow(P.x - position.x, 2) + pow(P.y - position.y, 2) + pow(P.z - position.z, 2));
	return min(1.0, 1.0/(constantTerm + linearTerm * d + quadraticTerm * pow(d, 2)));
}

glm::dvec3 PointLight::getColor() const
{
	printf("color\n");
	return color;
}

glm::dvec3 PointLight::getDirection(const glm::dvec3& P) const
{
	return glm::normalize(position - P);
}


glm::dvec3 PointLight::shadowAttenuation(const ray& r, const glm::dvec3& p) const
{
	printf("pointlight shadow atten\n");
	// YOUR CODE HERE:
	// You should implement shadow-handling code here.
	glm::dvec3 attenuation = r.getAtten();
	glm::dvec3 rayDirection = r.getDirection();
	glm::dvec3 directionOfLightFromRayIntersect = getDirection(p);
	ray shadow(p + RAY_EPSILON * (-rayDirection), directionOfLightFromRayIntersect, glm::dvec3(1,1,1), ray::SHADOW);
	isect shadowIntersect;
	bool intersect = scene->intersect(shadow, shadowIntersect);
	glm::dvec3 shadowIntersectPosition = shadow.at(shadowIntersect);
	glm::dvec3 directionOfLightFromShadowIntersect = getDirection(shadowIntersectPosition);
	bool in_front_of_light = glm::dot(directionOfLightFromShadowIntersect, directionOfLightFromRayIntersect) > 0;
	while (intersect && in_front_of_light) {
		Material m = shadowIntersect.getMaterial();
		if (m.Trans()) {
			bool leaving = glm::dot(directionOfLightFromShadowIntersect, shadowIntersect.getN()) > 0;
			if (leaving) {
				double d = glm::distance(shadow.getPosition(), shadowIntersectPosition);
				attenuation *= glm::pow(m.kt(shadowIntersect), glm::dvec3(d));
			}
		}
		else {
			attenuation = glm::dvec3(0.0, 0.0, 0.0);
		}
		in_front_of_light = glm::dot(getDirection(shadowIntersectPosition), directionOfLightFromShadowIntersect) > 0;
		directionOfLightFromShadowIntersect = getDirection(shadowIntersectPosition);
		shadow = ray(shadowIntersectPosition + RAY_EPSILON * directionOfLightFromShadowIntersect, directionOfLightFromShadowIntersect, glm::dvec3(1,1,1), ray::SHADOW);
		intersect = scene->intersect(shadow, shadowIntersect);
		shadowIntersectPosition = shadow.at(shadowIntersect);
	}
	printf("pointlight shadow atten finish\n");
	return attenuation;
}

#define VERBOSE 0

