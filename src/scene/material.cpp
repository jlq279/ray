#include "material.h"
#include "../ui/TraceUI.h"
#include "light.h"
#include "ray.h"
extern TraceUI* traceUI;

#include <glm/gtx/io.hpp>
#include <iostream>
#include "../fileio/images.h"

using namespace std;
extern bool debugMode;

Material::~Material()
{
}

// Apply the phong model to this point on the surface of the object, returning
// the color of that point.
glm::dvec3 Material::shade(Scene* scene, const ray& r, const isect& i) const
{
	// YOUR CODE HERE

	// For now, this method just returns the diffuse color of the object.
	// This gives a single matte color for every distinct surface in the
	// scene, and that's it.  Simple, but enough to get you started.
	// (It's also inconsistent with the phong model...)

	// Your mission is to fill in this method with the rest of the phong
	// shading model, including the contributions of all the light sources.
	// You will need to call both distanceAttenuation() and
	// shadowAttenuation()
	// somewhere in your code in order to compute shadows and light falloff.
	//	if( debugMode )
	//		std::cout << "Debugging Phong code..." << std::endl;

	// When you're iterating through the lights,
	// you'll want to use code that looks something
	// like this:
	//
	printf("scene %p\n", scene);
	glm::dvec3 rayDirection = r.getDirection();
	glm::dvec3 rayIntersectPosition = r.at(i);
	glm::dvec3 normal = i.getN();
	// glm::dvec3 emissive = ke(i);
	// glm::dvec3 ambient = ka(i) * scene->ambient();
	glm::dvec3 diffuse(0.0, 0.0, 0.0);
	// glm::dvec3 specular(0.0, 0.0, 0.0);
	// glm::dvec3 light(0.0, 0.0, 0.0);
	printf("about to iterate through lights\n");
	printf("rayIntersectPosition %f, %f, %f\n", rayIntersectPosition.x, rayIntersectPosition.y, rayIntersectPosition.z);

	for ( const auto& pLight : scene->getAllLights() ) {
		printf("plight %p\n", pLight.get());
		glm::dvec3 directionOfLightFromRayIntersect = pLight->getDirection(rayIntersectPosition);
		printf("directionOfLightFromRayIntersect %f, %f, %f\n", directionOfLightFromRayIntersect.x, directionOfLightFromRayIntersect.y, directionOfLightFromRayIntersect.z);
		printf("pLight->distanceAttenuation %f\n", pLight->distanceAttenuation(rayIntersectPosition));
		printf("shadow attenuation time, r %p\n", &r);
		glm::vec3 shadowAtten = pLight->shadowAttenuation(r, rayIntersectPosition);
		printf("pLight->shadowAttenuation %p\n", &shadowAtten);
		printf("pLight->getColor %f, %f, %f\n", pLight->getColor().r, pLight->getColor().g, pLight->getColor().b);
		auto iIn = pLight->distanceAttenuation(rayIntersectPosition) * pLight->shadowAttenuation(r, rayIntersectPosition) * pLight->getColor();
		printf("iIn %f, %f, %f\n", iIn.r, iIn.g, iIn.b);
		auto dot = glm::dot(directionOfLightFromRayIntersect, normal);
		if (Trans()) {
			dot = abs(dot);
		}
		else {
			dot = max(dot, 0.0);
		}
		glm::dvec3 pDiffuse = kd(i) * dot * iIn;
		diffuse += pDiffuse;
		// glm::dvec3 reflect = glm::reflect(-directionOfLightFromRayIntersect, normal);
		// auto maxDot = max(glm::dot(reflect, -rayDirection), 0.0);
		// glm::dvec3 pSpecular = ks(i) * pow(maxDot, shininess(i)) * iIn;
		// specular += pSpecular;
	}
	// return emissive + ambient + diffuse + specular;
	return diffuse;
}

TextureMap::TextureMap(string filename)
{
	data = readImage(filename.c_str(), width, height);
	if (data.empty()) {
		width = 0;
		height = 0;
		string error("Unable to load texture map '");
		error.append(filename);
		error.append("'.");
		throw TextureMapException(error);
	}
}

glm::dvec3 TextureMap::getMappedValue(const glm::dvec2& coord) const
{
	// YOUR CODE HERE
	//
	// In order to add texture mapping support to the
	// raytracer, you need to implement this function.
	// What this function should do is convert from
	// parametric space which is the unit square
	// [0, 1] x [0, 1] in 2-space to bitmap coordinates,
	// and use these to perform bilinear interpolation
	// of the values.
	double u = coord.x * (getWidth() - 1);
	double v = coord.y * (getHeight() - 1);
	double u1 = floor(u);
	double v1 = floor(v);
	double u2 = u1 + 1;
	double v2 = v1 + 1;
	double alpha = (u2 - u)/(u2 - u1);
	double beta = (u - u1)/(u2 - u1);
	glm::dvec3 a = getPixelAt(u1, v1);
	glm::dvec3 b = getPixelAt(u2, v1);
	glm::dvec3 c = getPixelAt(u2, v2);
	glm::dvec3 d = getPixelAt(u1, v2);
	return ((v2 - v)/(v2 - v1)) * (alpha * a + beta * b) + ((v - v1)/(v2 - v1)) * (alpha * d + beta * c);
}

glm::dvec3 TextureMap::getPixelAt(int x, int y) const
{
	// YOUR CODE HERE
	//
	// In order to add texture mapping support to the
	// raytracer, you need to implement this function.
	if (x < 0 || x >= getWidth() || y < 0 || y >= getHeight()) {
		return glm::dvec3(0, 0, 0);
	}
	auto index = 3 * (y * getWidth() + x);
	auto r = data.at(index);
	auto g = data.at(index + 1);
	auto b = data.at(index + 2);
	return glm::dvec3(r/255.0, g/255.0, b/255.0);
}

glm::dvec3 MaterialParameter::value(const isect& is) const
{
	if (0 != _textureMap)
		return _textureMap->getMappedValue(is.getUVCoordinates());
	else
		return _value;
}

double MaterialParameter::intensityValue(const isect& is) const
{
	if (0 != _textureMap) {
		glm::dvec3 value(
		        _textureMap->getMappedValue(is.getUVCoordinates()));
		return (0.299 * value[0]) + (0.587 * value[1]) +
		       (0.114 * value[2]);
	} else
		return (0.299 * _value[0]) + (0.587 * _value[1]) +
		       (0.114 * _value[2]);
}
