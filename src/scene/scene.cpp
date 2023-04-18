#include <cmath>

#include "scene.h"
#include "light.h"
#include "kdTree.h"
#include "../ui/TraceUI.h"
#include <glm/gtx/extended_min_max.hpp>
#include <iostream>
#include <glm/gtx/io.hpp>

using namespace std;

bool Geometry::intersect(ray& r, isect& i) const {
	// printf("geometry intersect\n");
	double tmin, tmax;
	if (hasBoundingBoxCapability() && !(bounds.intersect(r, tmin, tmax))) return false;
	// Transform the ray into the object's local coordinate space
	glm::dvec3 pos = transform->globalToLocalCoords(r.getPosition());
	glm::dvec3 dir = transform->globalToLocalCoords(r.getPosition() + r.getDirection()) - pos;
	double length = glm::length(dir);
	dir = glm::normalize(dir);
	// Backup World pos/dir, and switch to local pos/dir
	glm::dvec3 Wpos = r.getPosition();
	glm::dvec3 Wdir = r.getDirection();
	r.setPosition(pos);
	r.setDirection(dir);
	bool rtrn = false;
	if (intersectLocal(r, i))
	{
		// Transform the intersection point & normal returned back into global space.
		i.setN(transform->localToGlobalCoordsNormal(i.getN()));
		i.setT(i.getT()/length);
		rtrn = true;
	}
	// Restore World pos/dir
	r.setPosition(Wpos);
	r.setDirection(Wdir);
	return rtrn;
}

bool Geometry::hasBoundingBoxCapability() const {
	// by default, primitives do not have to specify a bounding box.
	// If this method returns true for a primitive, then either the ComputeBoundingBox() or
    // the ComputeLocalBoundingBox() method must be implemented.

	// If no bounding box capability is supported for an object, that object will
	// be checked against every single ray drawn.  This should be avoided whenever possible,
	// but this possibility exists so that new primitives will not have to have bounding
	// boxes implemented for them.
	return false;
}

void Geometry::ComputeBoundingBox() {
    // take the object's local bounding box, transform all 8 points on it,
    // and use those to find a new bounding box.

    BoundingBox localBounds = ComputeLocalBoundingBox();
        
    glm::dvec3 min = localBounds.getMin();
    glm::dvec3 max = localBounds.getMax();

    glm::dvec4 v, newMax, newMin;

    v = transform->localToGlobalCoords( glm::dvec4(min[0], min[1], min[2], 1) );
    newMax = v;
    newMin = v;
    v = transform->localToGlobalCoords( glm::dvec4(max[0], min[1], min[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(min[0], max[1], min[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(max[0], max[1], min[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(min[0], min[1], max[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(max[0], min[1], max[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(min[0], max[1], max[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(max[0], max[1], max[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
		
    bounds.setMax(glm::dvec3(newMax));
    bounds.setMin(glm::dvec3(newMin));
}

Scene::Scene()
{
	ambientIntensity = glm::dvec3(0, 0, 0);
}

Scene::~Scene()
{
}

void Scene::add(Geometry* obj) {
	obj->ComputeBoundingBox();
	sceneBounds.merge(obj->getBoundingBox());
	objects.emplace_back(obj);
}

void Scene::add(Light* light)
{
	lights.emplace_back(light);
}

bool traverse(KdTree<Geometry>* tree, ray& r, isect& i) {
	double tMin;
	double tMax;
	BoundingBox bound = tree->getBound();
	if (bound.intersect(r, tMin, tMax)) {
		KdTree<Geometry>* left = tree->getLeft();
		KdTree<Geometry>* right = tree->getRight();
		if (tree->isLeaf()) {
			bool have_one = false; 
			for (auto& object : tree->getObjects()) {
				auto b = object.get()->getBoundingBox();
				isect cur;
				if( object->intersect(r, cur) ) {
					if((!have_one && i.getT() == 0) || (cur.getT() < i.getT())) {
						i = cur;
						have_one = true;
					}
				}
			}
			return have_one;
		}
		bool leftIntersect = left == NULL ? false : traverse(left, r, i);
		bool rightIntersect = right == NULL ? false : traverse(right, r, i);
		return leftIntersect || rightIntersect;
	}
	else {
		return false;
	}
}

// Get any intersection with an object.  Return information about the 
// intersection through the reference parameter.
bool Scene::intersect(ray& r, isect& i) const {
	bool have_one = traverse(kdTree, r, i);
	if(!have_one)
		i.setT(1000.0);
	// if debugging,
	if (TraceUI::m_debug)
	{
		addToIntersectCache(std::make_pair(new ray(r), new isect(i)));
	}
	return have_one;
}

TextureMap* Scene::getTexture(string name) {
	auto itr = textureCache.find(name);
	if (itr == textureCache.end()) {
		textureCache[name].reset(new TextureMap(name));
		return textureCache[name].get();
	}
	return itr->second.get();
}

bool compareX(std::unique_ptr<Geometry>& a, std::unique_ptr<Geometry>& b) {
	if (a->getBoundingBox().getMin().x == b->getBoundingBox().getMin().x) {
		if (a->getBoundingBox().getMax().x == b->getBoundingBox().getMax().x) {
			return sqrt(pow(a->getBoundingBox().getMax().x, 2) + pow(a->getBoundingBox().getMax().y, 2) + pow(a->getBoundingBox().getMax().z, 2)) < sqrt(pow(b->getBoundingBox().getMax().x, 2) + pow(b->getBoundingBox().getMax().y, 2) + pow(b->getBoundingBox().getMax().z, 2));
		}
		return a->getBoundingBox().getMax().x < b->getBoundingBox().getMax().x;
	}
	return a->getBoundingBox().getMin().x < b->getBoundingBox().getMin().x;
}

bool compareY(std::unique_ptr<Geometry>& a, std::unique_ptr<Geometry>& b) {
	if (a->getBoundingBox().getMin().y == b->getBoundingBox().getMin().y) {
		if (a->getBoundingBox().getMax().y < b->getBoundingBox().getMax().y) {
			return sqrt(pow(a->getBoundingBox().getMax().x, 2) + pow(a->getBoundingBox().getMax().y, 2) + pow(a->getBoundingBox().getMax().z, 2)) < sqrt(pow(b->getBoundingBox().getMax().x, 2) + pow(b->getBoundingBox().getMax().y, 2) + pow(b->getBoundingBox().getMax().z, 2));
		}
		return a->getBoundingBox().getMax().y < b->getBoundingBox().getMax().y;
	}
	return a->getBoundingBox().getMin().y < b->getBoundingBox().getMin().y;
}

bool compareZ(std::unique_ptr<Geometry>& a, std::unique_ptr<Geometry>& b) {
	if (a->getBoundingBox().getMin().z == b->getBoundingBox().getMin().z) {
		if (a->getBoundingBox().getMax().z == b->getBoundingBox().getMax().z) {
			return sqrt(pow(a->getBoundingBox().getMax().x, 2) + pow(a->getBoundingBox().getMax().y, 2) + pow(a->getBoundingBox().getMax().z, 2)) < sqrt(pow(b->getBoundingBox().getMax().x, 2) + pow(b->getBoundingBox().getMax().y, 2) + pow(b->getBoundingBox().getMax().z, 2));
		}
		return a->getBoundingBox().getMax().z < b->getBoundingBox().getMax().z;
	}
	return a->getBoundingBox().getMin().z < b->getBoundingBox().getMin().z;
}

void buildKdTreeR(KdTree<Geometry>* tree, int depth, int leafSize) {
	int numObjects = tree->getNumObjects();
	BoundingBox bound = tree->getBound();
	auto min = bound.getMin();
	auto max = bound.getMax();
	auto minX = min.x;
	auto maxX = max.x;
	auto minY = min.y;
	auto maxY = max.y;
	auto minZ = min.z;
	auto maxZ = max.z;
	auto difference = max - min;
	if ((depth == 0 || numObjects <= leafSize)) {
		return;
	}
	KdTree<Geometry>* left;
	KdTree<Geometry>* right;
	if (difference.x >= difference.y && difference.x >= difference.z) {
		sort(tree->getObjects().begin(), tree->getObjects().end(), compareX);
	}
	else if (difference.y >= difference.x && difference.y >= difference.z) {
		sort(tree->getObjects().begin(), tree->getObjects().end(), compareY);
	}
	else {
		sort(tree->getObjects().begin(), tree->getObjects().end(), compareZ);
	}
	auto half = tree->getObjects().size()/2;
	left = new KdTree<Geometry>();
	std::vector<std::unique_ptr<Geometry>> leftObjects;
	right = new KdTree<Geometry>();
	std::vector<std::unique_ptr<Geometry>> rightObjects;
	for (auto object = tree->getObjects().begin(); object < tree->getObjects().end(); object++) {
		if (object < tree->getObjects().begin() + half) {
			left->addObject(object->get());
			left->updateBound(object->get()->getBoundingBox());
		}
		else {
			right->addObject(object->get());
			right->updateBound(object->get()->getBoundingBox());
		}
	}
	if (left->getNumObjects() > 0) {
		tree->setLeft(left);
		buildKdTreeR(tree->getLeft(), depth - 1, leafSize);
	}
	if (right->getNumObjects() > 0) {
		tree->setRight(right);
		buildKdTreeR(tree->getRight(), depth - 1, leafSize);
	}
}

void Scene::buildKdTree(int depth, int leafSize) {
	kdTree = new KdTree<Geometry>(sceneBounds.getMin(), sceneBounds.getMax());
	for ( auto& object : objects ) {
		auto bb = object.get()->getBoundingBox();
		kdTree->addObject(object.get());
		kdTree->updateBound(bb.getMin(), bb.getMax());
	}
	buildKdTreeR(kdTree, depth, leafSize);
}
