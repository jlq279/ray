#ifndef __RAYTRACER_H__
#define __RAYTRACER_H__

#define MAX_THREADS 32

// The main ray tracer.

#include <time.h>
#include <glm/vec3.hpp>
#include <queue>
#include <thread>
#include "scene/cubeMap.h"
#include "scene/ray.h"
#include <mutex>

class Scene;
class Pixel {
public:
	Pixel(int i, int j, unsigned char* ptr) : ix(i), jy(j), value(ptr) {}

	int ix;
	int jy;
	unsigned char* value;
};


class RayTracer {
public:
	RayTracer();
	~RayTracer();

	glm::dvec3 tracePixel(int i, int j);
	glm::dvec3 traceRay(ray& r, const glm::dvec3& thresh, int depth,
	                    double& length);
	glm::dvec3 traceOffset(double x, double y, double offset);

	glm::dvec3 getPixel(int i, int j);
	void setPixel(int i, int j, glm::dvec3 color);
	void getBuffer(unsigned char*& buf, int& w, int& h);
	double aspectRatio();

	void launchThread();
	void traceImage(int w, int h);
	int aaImage();
	int assAAImage();
	glm::dvec3 adaptiveSupersampling(double x, double y, double dimension, int depth);
	bool checkRender();
	void waitRender();

	void traceSetup(int w, int h);

	bool loadScene(const char* fn);
	bool sceneLoaded() { return scene != 0; }

	void setReady(bool ready) { m_bBufferReady = ready; }
	bool isReady() const { return m_bBufferReady; }

	const Scene& getScene() { return *scene; }

	bool stopTrace;

private:
	glm::dvec3 trace(double x, double y);

	std::vector<unsigned char> buffer;
	int buffer_width, buffer_height;
	int bufferSize;
	unsigned int threads;
	int block_size;
	double thresh;
	double aaThresh;
	int samples;
	int rays;
	double m_3dOffset; 
	std::unique_ptr<Scene> scene;
	std::vector<std::thread> threads_arr;
	std::vector<int> threads_indices;
	std::vector<bool> threads_finished;
	int w;
	int h;
	std::mutex myMutex;


	bool m_bBufferReady;
	bool renderReady;
};

#endif // __RAYTRACER_H__
