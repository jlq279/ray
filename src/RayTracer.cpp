// The main ray tracer.

#pragma warning (disable: 4786)

#include "RayTracer.h"
#include "scene/light.h"
#include "scene/material.h"
#include "scene/ray.h"

#include "parser/Tokenizer.h"
#include "parser/Parser.h"

#include "ui/TraceUI.h"
#include <cmath>
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>
#include <string.h> // for memset
#include <random>

#include <iostream>
#include <fstream>

using namespace std;
extern TraceUI* traceUI;

// Use this variable to decide if you want to print out
// debugging messages.  Gets set in the "trace single ray" mode
// in TraceGLWindow, for example.
bool debugMode = false;

// Trace a top-level ray through pixel(i,j), i.e. normalized window coordinates (x,y),
// through the projection plane, and out into the scene.  All we do is
// enter the main ray-tracing method, getting things started by plugging
// in an initial ray weight of (0.0,0.0,0.0) and an initial recursion depth of 0.

glm::dvec3 RayTracer::trace(double x, double y)
{
	// Clear out the ray cache in the scene for debugging purposes,
	if (TraceUI::m_debug)
	{
		scene->clearIntersectCache();		
	}

	ray r(glm::dvec3(0,0,0), glm::dvec3(0,0,0), glm::dvec3(1,1,1), ray::VISIBILITY);
	scene->getCamera().rayThrough(x,y,r);
	double dummy;
	glm::dvec3 ret = traceRay(r, glm::dvec3(1.0,1.0,1.0), traceUI->getDepth(), dummy);
	ret = glm::clamp(ret, 0.0, 1.0);
	return ret;
}

glm::dvec3 RayTracer::tracePixel(int i, int j)
{
	glm::dvec3 col(0,0,0);

	if( ! sceneLoaded() ) return col;

	double x = double(i)/double(buffer_width);
	double y = double(j)/double(buffer_height);

	unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;
	col = trace(x, y);

	pixel[0] = (int)( 255.0 * col[0]);
	pixel[1] = (int)( 255.0 * col[1]);
	pixel[2] = (int)( 255.0 * col[2]);
	return col;
}

#define VERBOSE 0

void createCoordinateSystem(glm::dvec3 &N, glm::dvec3 &Nt, glm::dvec3 &Nb) 
{ 
	if (std::fabs(N.x) > std::fabs(N.y))
	{
		float scale = sqrtf(N.x * N.x + N.z * N.z);
		Nt = glm::dvec3(N.z/scale, 0, -N.x/scale);
	}
	else
	{
		float scale = sqrtf(N.y * N.y + N.z * N.z);
		Nt = glm::dvec3(0, -N.z/scale, N.y/scale);
	}
	Nb = glm::cross(N, Nt);
}

glm::dvec3 uniformSampleHemisphere(const float &r1, const float &r2)
{
    // cos(theta) = r1 = y
    // cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
    float sinTheta = sqrtf(1 - r1 * r1);
    float phi = 2 * M_PI * r2;
    float x = sinTheta * cosf(phi);
    float z = sinTheta * sinf(phi);
    return glm::dvec3(x, r1, z);
}

// Do recursive ray tracing!  You'll want to insert a lot of code here
// (or places called from here) to handle reflection, refraction, etc etc.
glm::dvec3 RayTracer::traceRay(ray& r, const glm::dvec3& thresh, int depth, double& t )
{
	isect i;
	glm::dvec3 colorC;
#if VERBOSE
	std::cerr << "== current depth: " << depth << std::endl;
#endif

	if(scene->intersect(r, i)) {
		// YOUR CODE HERE

		// An intersection occurred!  We've got work to do.  For now,
		// this code gets the material for the surface that was intersected,
		// and asks that material to provide a color for the ray.

		// This is a great place to insert code for recursive ray tracing.
		// Instead of just returning the result of shade(), add some
		// more steps: add in the contributions from reflected and refracted
		// rays.

		// printf("traceray intersect\n");
		
		const Material& m = i.getMaterial();
		// printf("material\n");
		glm::dvec3 direct = m.shade(scene.get(), r, i);
		// printf("colorC %f, %f, %f\n", colorC.r, colorC.g, colorC.b);
		glm::dvec3 intersectPosition = r.at(i);
		// printf("intersectPosition %f, %f, %f\n", intersectPosition.x, intersectPosition.y, intersectPosition.z);
		glm::dvec3 n = i.getN();
		// printf("n %f, %f, %f\n", n.x, n.y, n.z);
		// if (m.Recur() && depth > 0) {
		// 	glm::dvec3 rayPosition = r.getPosition();
		// 	glm::dvec3 rayDirection = r.getDirection();
		// 	bool leaving = glm::dot(rayDirection, n) > 0;
		// 	if (m.Refl()) {
		// 		if (leaving) {
		// 			glm::dvec3 reflect = glm::reflect(rayDirection, -n);
		// 			ray reflectedRay(intersectPosition + RAY_EPSILON * (-rayDirection), reflect, glm::dvec3(1,1,1), ray::REFLECTION);
		// 			glm::dvec3 reflectedColor = traceRay(reflectedRay, thresh, depth - 1, t);
		// 			direct += reflectedColor * m.kr(i);
		// 		}
		// 		else {
		// 			glm::dvec3 reflect = glm::reflect(rayDirection, n);
		// 			ray reflectedRay(intersectPosition + RAY_EPSILON * (-rayDirection), reflect, glm::dvec3(1,1,1), ray::REFLECTION);
		// 			glm::dvec3 reflectedColor = traceRay(reflectedRay, thresh, depth - 1, t);
		// 			direct += reflectedColor * m.kr(i);
		// 		}
				
		// 	}
		// 	if (m.Trans()) {
		// 		if (leaving) {
		// 			double eta = m.index(i);
		// 			double d = glm::distance(rayPosition, intersectPosition);
		// 			glm::dvec3 refract = glm::refract(rayDirection, -n, eta);
		// 			if (refract.x == 0 && refract.y == 0 && refract.z == 0) {
		// 				glm::dvec3 reflect = glm::reflect(rayDirection, -n);
		// 				ray reflectedRay(intersectPosition + RAY_EPSILON * (-rayDirection), reflect, glm::dvec3(1,1,1), ray::REFLECTION);
		// 				glm::dvec3 reflectedColor = traceRay(reflectedRay, thresh, depth - 1, t);
		// 				direct += reflectedColor * m.kr(i);
		// 			}
		// 			else {
		// 				ray refractedRay(intersectPosition + RAY_EPSILON * rayDirection, refract, glm::dvec3(1,1,1), ray::REFRACTION);
		// 				glm::dvec3 refractedColor = traceRay(refractedRay, thresh, depth - 1, t);
		// 				direct += refractedColor * glm::pow(m.kt(i), glm::dvec3(d));
		// 			}
		// 		}
		// 		else {
		// 			double eta = 1.0/m.index(i);
		// 			glm::dvec3 refract = glm::refract(rayDirection, n, eta);
		// 			if (refract.x == 0 && refract.y == 0 && refract.z == 0) {
		// 				glm::dvec3 reflect = glm::reflect(rayDirection, n);
		// 				ray reflectedRay(intersectPosition + RAY_EPSILON * (-rayDirection), reflect, glm::dvec3(1,1,1), ray::REFLECTION);
		// 				glm::dvec3 reflectedColor = traceRay(reflectedRay, thresh, depth - 1, t);
		// 				direct += reflectedColor * m.kr(i);
		// 			}
		// 			else {
		// 				ray refractedRay(intersectPosition + RAY_EPSILON * rayDirection, refract, glm::dvec3(1,1,1), ray::REFRACTION);
		// 				glm::dvec3 refractedColor = traceRay(refractedRay, thresh, depth - 1, t);
		// 				direct += refractedColor;
		// 			}
		// 		}
				
		// 	}
			
		// }
		if (depth > 0)
		{
			glm::dvec3 nt(0.0, 0.0, 0.0);
			glm::dvec3 nb(0.0, 0.0, 0.0);
			createCoordinateSystem(n, nt, nb);
			// printf("nt %f, %f, %f\n", nt.x, nt.y, nt.z);
			// printf("nb %f, %f, %f\n", nb.x, nb.y, nb.z);
			float pdf = 1 / (2 * M_PI);
			//std::random_dev dev;
			//std::mt19937 rng(dev());

			std::default_random_engine generator;
			std::uniform_real_distribution<float> distribution(0, 1);
			uint32_t N = 16;
			glm::dvec3 indirect(0.0, 0.0, 0.0);
			for (uint32_t i = 0; i < N; ++i) {
				float r1 = ((double) rand() / RAND_MAX);
				float r2 = ((double) rand() / RAND_MAX);
				// printf("r1 %f, r2 %f\n", r1, r2);
				glm::dvec3 sample = uniformSampleHemisphere(r1, r2);
				// printf("sample %f, %f, %f\n", sample.x, sample.y, sample.z);
				glm::dvec3 sampleWorld( 
					sample.x * nb.x + sample.y * n.x + sample.z * nt.x,
					sample.x * nb.y + sample.y * n.y + sample.z * nt.y,
					sample.x * nb.z + sample.y * n.z + sample.z * nt.z);
				// printf("sampleWorld %f, %f, %f\n", sampleWorld.x, sampleWorld.y, sampleWorld.z);
				ray sampleRay(intersectPosition + sampleWorld * RAY_EPSILON, sampleWorld, glm::dvec3(1,1,1), ray::REFLECTION);
				// printf("sampleRay %f, %f, %f\n", sampleRay.getPosition().x, sampleRay.getPosition().y, sampleRay.getPosition().z);
				glm::dvec3 traced = traceRay(sampleRay, thresh, depth - 1, t);
				// printf("traced %f, %f, %f\n", traced.r, traced.g, traced.b);
				indirect += glm::dvec3(r1 * traced.r, r1 * traced.g, r1 * traced.b);
			}
			indirect /= (float) N;
			colorC = (direct + glm::dvec3(2 * indirect.r, 2 * indirect.g, 2 * indirect.b)) * m.kd(i);
		}
		else {
			colorC = direct * m.kd(i);
		}
	} else {
		// No intersection.  This ray travels to infinity, so we color
		// it according to the background color, which in this (simple) case
		// is just black.
		//
		// FIXME: Add CubeMap support here.
		// TIPS: CubeMap object can be fetched from traceUI->getCubeMap();
		//       Check traceUI->cubeMap() to see if cubeMap is loaded
		//       and enabled.
		// printf("cube\n");
		if (traceUI->cubeMap()) {
			auto cubeMap = traceUI->getCubeMap();
			colorC = cubeMap->getColor(r);
		}
		else {
			colorC = glm::dvec3(0.0, 0.0, 0.0);
		}
	}
#if VERBOSE
	std::cerr << "== depth: " << depth+1 << " done, returning: " << colorC << std::endl;
#endif
	return colorC;
}

RayTracer::RayTracer()
	: scene(nullptr), buffer(0), thresh(0), buffer_width(0), buffer_height(0), m_bBufferReady(false)
{
}

RayTracer::~RayTracer()
{
}

void RayTracer::getBuffer( unsigned char *&buf, int &w, int &h )
{
	buf = buffer.data();
	w = buffer_width;
	h = buffer_height;
}

double RayTracer::aspectRatio()
{
	return sceneLoaded() ? scene->getCamera().getAspectRatio() : 1;
}

bool RayTracer::loadScene(const char* fn)
{
	ifstream ifs(fn);
	if( !ifs ) {
		string msg( "Error: couldn't read scene file " );
		msg.append( fn );
		traceUI->alert( msg );
		return false;
	}

	// Strip off filename, leaving only the path:
	string path( fn );
	if (path.find_last_of( "\\/" ) == string::npos)
		path = ".";
	else
		path = path.substr(0, path.find_last_of( "\\/" ));

	// Call this with 'true' for debug output from the tokenizer
	Tokenizer tokenizer( ifs, false );
	Parser parser( tokenizer, path );
	try {
		scene.reset(parser.parseScene());
	}
	catch( SyntaxErrorException& pe ) {
		traceUI->alert( pe.formattedMessage() );
		return false;
	} catch( ParserException& pe ) {
		string msg( "Parser: fatal exception " );
		msg.append( pe.message() );
		traceUI->alert( msg );
		return false;
	} catch( TextureMapException e ) {
		string msg( "Texture mapping exception: " );
		msg.append( e.message() );
		traceUI->alert( msg );
		return false;
	}

	if (!sceneLoaded())
		return false;

	return true;
}

void RayTracer::traceSetup(int w, int h)
{
	size_t newBufferSize = w * h * 3;
	if (newBufferSize != buffer.size()) {
		bufferSize = newBufferSize;
		buffer.resize(bufferSize);
	}
	buffer_width = w;
	buffer_height = h;
	std::fill(buffer.begin(), buffer.end(), 0);
	m_bBufferReady = true;

	/*
	 * Sync with TraceUI
	 */

	threads = traceUI->getThreads();
	block_size = traceUI->getBlockSize();
	thresh = traceUI->getThreshold();
	samples = traceUI->getSuperSamples();
	aaThresh = traceUI->getAaThreshold();
	m_3dOffset = traceUI->get3dOffset();

	// YOUR CODE HERE
	// FIXME: Additional initializations

	if (sceneLoaded()) {
		auto depth = traceUI->getMaxDepth();
		auto leafSize = traceUI->getLeafSize();
		scene->buildKdTree(depth, leafSize);
	}
}

glm::dvec3 RayTracer::traceOffset(double x, double y, double offset)
{
	ray r(glm::dvec3(0,0,0), glm::dvec3(0,0,0), glm::dvec3(1,1,1), ray::VISIBILITY);
	glm::dvec3 up = scene->getCamera().up;
	glm::dvec3 cross = glm::cross(scene->getCamera().getLook(), up);
	glm::dvec3 horizontal = glm::normalize(cross);
	const glm::dvec3 eye(scene->getCamera().getEye());
	glm::dvec3 leftEye(eye.x, eye.y, eye.z);
	leftEye = leftEye - offset * horizontal;
 	scene->getCamera().setEye(leftEye);
	glm::dvec3 left = trace(x, y);
	glm::dvec3 rightEye(eye.x, eye.y, eye.z);
	rightEye = rightEye + offset * horizontal;
	scene->getCamera().setEye(rightEye);
	glm::dvec3 right = trace(x, y);
	scene->getCamera().setEye(eye);
	glm::dvec3 filter = glm::dvec3(0.299, 0.587, 0.114);
	glm::dvec3 red = left * filter;
	glm::dvec3 green = right * filter;
	glm::dvec3 blue = right * filter;
	return glm::dvec3(red.r + red.g + red.b, green.r + green.g + green.b, blue.r + blue.g + blue.b);
}

/*
 * RayTracer::traceImage
 *
 *	Trace the image and store the pixel data in RayTracer::buffer.
 *
 *	Arguments:
 *		w:	width of the image buffer
 *		h:	height of the image buffer
 *
 */
void RayTracer::traceImage(int w, int h)
{
	// Always call traceSetup before rendering anything.
	traceSetup(w,h);

	// YOUR CODE HERE
	// FIXME: Start one or more threads for ray tracing
	//
	// TIPS: Ideally, the traceImage should be executed asynchronously,
	//       i.e. returns IMMEDIATELY after working threads are launched.
	//
	//       An asynchronous traceImage lets the GUI update your results
	//       while rendering.
	if (traceUI->m_3dSwitch()) {
		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				glm::dvec3 color = traceOffset(double(i)/double(buffer_width), double(j)/double(buffer_height), m_3dOffset/100);
				unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;
				pixel[0] = (int)( 255.0 * color[0]);
				pixel[1] = (int)( 255.0 * color[1]);
				pixel[2] = (int)( 255.0 * color[2]);
			}
		}
	}
	else {
		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				tracePixel(i, j);
			}
		}
	}
}

glm::dvec3 RayTracer::adaptiveSupersampling(double x, double y, double dimension, int depth) {
	glm::dvec3 center = trace(x, y);
	glm::dvec3 top_left = trace(x - dimension/2, y);
	glm::dvec3 top_right = trace(x + dimension/2, y);
	glm::dvec3 bottom_left = trace(x - dimension/2, y + dimension/2);
	glm::dvec3 bottom_right = trace(x + dimension/2, y + dimension/2);
	// rays += 5;
	glm::dvec3 average = (center + top_left + top_right + bottom_left + bottom_right)/5.0;
	if (depth == 0) {
		return average;
	}
	glm::dvec3 d0 = glm::abs(center - average);
	glm::dvec3 d1 = glm::abs(top_left - average);
	glm::dvec3 d2 = glm::abs(top_right - average);
	glm::dvec3 d3 = glm::abs(bottom_left - average);
	glm::dvec3 d4 = glm::abs(bottom_right - average);

	if ((d0.x < aaThresh && d0.y < aaThresh && d0.z < aaThresh) &&
		(d1.x < aaThresh && d1.y < aaThresh && d1.z < aaThresh) &&
		(d2.x < aaThresh && d2.y < aaThresh && d2.z < aaThresh) &&
		(d3.x < aaThresh && d3.y < aaThresh && d3.z < aaThresh) &&
		(d4.x < aaThresh && d4.y < aaThresh && d4.z < aaThresh)) {
			return average;
	}
	else {
		double subdimension = dimension/samples;
		glm::dvec3 subaverage(0.0, 0.0, 0.0);
		for (double k = 0; k < samples; k++) {
			for (double l = 0; l < samples; l++) {
				double x1 = x - dimension/2 + k * subdimension + subdimension/2;
				double y1 = y - dimension/2 + l * subdimension + subdimension/2;
				if (glm::abs(x1 - x) < 0.0001 && glm::abs(y1 - y) < 0.0001) {
					return average;
				}
				subaverage += adaptiveSupersampling(x1, y1, subdimension, depth - 1);
			}
		}
		return subaverage/pow(samples, 2);
	}
}

int RayTracer::assAAImage() {
	for (int i = 0; i < buffer_width; i++) {
		for (int j = 0; j < buffer_height; j++) {
			glm::dvec3 average(0.0, 0.0, 0.0);
			double pixel_dimension = 1.0/buffer_width;
			double x = i * pixel_dimension + pixel_dimension/2;
			double y = j * pixel_dimension + pixel_dimension/2;
			// rays = 0;
			glm::dvec3 color = adaptiveSupersampling(x, y, pixel_dimension, 10);
			unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;
			pixel[0] = (int)( 255.0 * color[0]);
			pixel[1] = (int)( 255.0 * color[1]);
			pixel[2] = (int)( 255.0 * color[2]);
			// pixel[0] = rays;
			// pixel[1] = rays;
			// pixel[2] = rays;
		}
	}
	return 0;
}

int RayTracer::aaImage()
{
	// YOUR CODE HERE
	// FIXME: Implement Anti-aliasing here
	//
	// TIP: samples and aaThresh have been synchronized with TraceUI by
	//      RayTracer::traceSetup() function
	if (traceUI->assAASwitch()) {
		assAAImage();
	}
	else {
		for (int i = 0; i < buffer_width; i++) {
			for (int j = 0; j < buffer_height; j++) {
				glm::dvec3 average(0.0, 0.0, 0.0);
				double pixel_dimension = 1.0/buffer_width;
				double subpixel_dimension = pixel_dimension/samples;
				for (double k = 0; k < samples; k++) {
					for (double l = 0; l < samples; l++) {
						double x = i * pixel_dimension + (k + 0.5) * subpixel_dimension;
						double y = j * pixel_dimension + (l + 0.5) * subpixel_dimension;
						glm::dvec3 traced = trace(x, y);
						average += traced;
					}
				}
				average /= pow(samples, 2);
				unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;
				pixel[0] = (int)( 255.0 * average[0]);
				pixel[1] = (int)( 255.0 * average[1]);
				pixel[2] = (int)( 255.0 * average[2]);

			}
		}
	}
	return 0;
}

bool RayTracer::checkRender()
{
	// YOUR CODE HERE
	// FIXME: Return true if tracing is done.
	//        This is a helper routine for GUI.
	//
	// TIPS: Introduce an array to track the status of each worker thread.
	//       This array is maintained by the worker threads.
	return true;
}

void RayTracer::waitRender()
{
	// YOUR CODE HERE
	// FIXME: Wait until the rendering process is done.
	//        This function is essential if you are using an asynchronous
	//        traceImage implementation.
	//
	// TIPS: Join all worker threads here.
}


glm::dvec3 RayTracer::getPixel(int i, int j)
{
	unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;
	return glm::dvec3((double)pixel[0]/255.0, (double)pixel[1]/255.0, (double)pixel[2]/255.0);
}

void RayTracer::setPixel(int i, int j, glm::dvec3 color)
{
	unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;

	pixel[0] = (int)( 255.0 * color[0]);
	pixel[1] = (int)( 255.0 * color[1]);
	pixel[2] = (int)( 255.0 * color[2]);
}

