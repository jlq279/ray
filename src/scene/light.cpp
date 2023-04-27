#include <cmath>
#include <iostream>

#include "light.h"
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>
#include <random>

using namespace std;

double DirectionalLight::distanceAttenuation(const glm::dvec3& P) const
{
	// distance to light is infinite, so f(di) goes to 0.  Return 1.
	// printf("test");
	return 1.0;
}


glm::dvec3 DirectionalLight::shadowAttenuation(const ray& r, const glm::dvec3& p) const
{
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
	return attenuation;
}

glm::dvec3 DirectionalLight::getColor() const
{
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
	double d = sqrt(pow(P.x - position.x, 2) + pow(P.y - position.y, 2) + pow(P.z - position.z, 2));
	return min(1.0, 1.0/(constantTerm + linearTerm * d + quadraticTerm * pow(d, 2)));
}

glm::dvec3 PointLight::getColor() const
{
	return color;
}

glm::dvec3 PointLight::getDirection(const glm::dvec3& P) const
{
	return glm::normalize(position - P);
}

glm::dvec3 sampleShadow()
{
    // cos(theta) = r1 = y
    // cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
	float r1 = ((double) rand() / RAND_MAX);
	float r2 = ((double) rand() / RAND_MAX);
    float sinTheta = sqrtf(1 - r1 * r1);
    float phi = 2 * M_PI * r2;
    float x = sinTheta * cosf(phi);
    float z = sinTheta * sinf(phi);
    return glm::dvec3(x * 0.5, r1 * 0.5, z * 0.5);
}

glm::dvec3 PointLight::shadowAttenuation(const ray& r, const glm::dvec3& p) const
{
	// YOUR CODE HERE:
	// You should implement shadow-handling code here.
	glm::dvec3 rayDirection = r.getDirection();
	glm::dvec3 attenuation = r.getAtten();
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
			break;
		}
		in_front_of_light = glm::dot(getDirection(shadowIntersectPosition), directionOfLightFromShadowIntersect) > 0;
		directionOfLightFromShadowIntersect = getDirection(shadowIntersectPosition);
		shadow = ray(shadowIntersectPosition + RAY_EPSILON * directionOfLightFromShadowIntersect, directionOfLightFromShadowIntersect, glm::dvec3(1,1,1), ray::SHADOW);
		intersect = scene->intersect(shadow, shadowIntersect);
		shadowIntersectPosition = shadow.at(shadowIntersect);
	}
	return attenuation;
}

double SpotLight::distanceAttenuation(const glm::dvec3& P) const
{

	// YOUR CODE HERE

	// You'll need to modify this method to attenuate the intensity 
	// of the light based on the distance between the source and the 
	// point P.  For now, we assume no attenuation and just return 1.0
	// f(d) = min( 1, 1/( a + b d + c d^2 ) )
	// float constantTerm;		// a
	// float linearTerm;		// b
	// float quadraticTerm;		// c
	double d = sqrt(pow(P.x - position.x, 2) + pow(P.y - position.y, 2) + pow(P.z - position.z, 2));
	return min(1.0, 1.0/(constantTerm + linearTerm * d + quadraticTerm * pow(d, 2)));
}

glm::dvec3 SpotLight::getColor() const
{
	return color;
}

glm::dvec3 SpotLight::getDirection(const glm::dvec3& P) const
{
	return glm::normalize(position - P);
}

glm::dvec3 samplePoint()
{
    // cos(theta) = r1 = y
    // cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
	float r1 = ((double) rand() / RAND_MAX);
	float r2 = ((double) rand() / RAND_MAX);
    float sinTheta = sqrtf(1 - r1 * r1);
    float phi = 2 * M_PI * r2;
    float x = sinTheta * cosf(phi);
    float z = sinTheta * sinf(phi);
    return glm::dvec3(x * 0.5, r1 * 0.5, z * 0.5);
}

glm::dvec3 SpotLight::shadowAttenuation(const ray& r, const glm::dvec3& p) const
{
	// YOUR CODE HERE:
	// You should implement shadow-handling code here.
	glm::dvec3 rayDirection = r.getDirection();
	uint32_t N = 16;
	glm::dvec3 shadowAttenuation(0.0, 0.0, 0.0);
	for (uint32_t i = 0; i < N; ++i) {
		glm::dvec3 attenuation = r.getAtten();
		glm::dvec3 lightPosition = position + sampleShadow();
		glm::dvec3 directionOfLightFromRayIntersect = glm::normalize(lightPosition - p);
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
				break;
			}
			in_front_of_light = glm::dot(getDirection(shadowIntersectPosition), directionOfLightFromShadowIntersect) > 0;
			directionOfLightFromShadowIntersect = getDirection(shadowIntersectPosition);
			shadow = ray(shadowIntersectPosition + RAY_EPSILON * directionOfLightFromShadowIntersect, directionOfLightFromShadowIntersect, glm::dvec3(1,1,1), ray::SHADOW);
			intersect = scene->intersect(shadow, shadowIntersect);
			shadowIntersectPosition = shadow.at(shadowIntersect);
		}
		shadowAttenuation += attenuation;
	}
	return glm::dvec3(shadowAttenuation.x/N, shadowAttenuation.y/N, shadowAttenuation.z/N);
}

#define VERBOSE 0

