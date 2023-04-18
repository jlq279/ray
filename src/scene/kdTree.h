#pragma once
#include "bbox.h"
#include <vector>
// #include <glm/gtc/matrix_transform.hpp>
// #include <glm/vec3.hpp>

// Note: you can put kd-tree here
template <typename Obj>
class KdTree {
	BoundingBox bound;
	std::vector<std::unique_ptr<Obj>> objects = {};
    KdTree* left = 0;
    KdTree* right = 0;
public:
    KdTree() {}
    KdTree(BoundingBox b) : bound(b) {}
    KdTree(glm::dvec3 bMin, glm::dvec3 bMax) : bound(BoundingBox(bMin, bMax)) {}
	BoundingBox getBound() { return bound; }
	KdTree* getLeft() { return left; }
	KdTree* getRight() { return right; }
	void setLeft(KdTree* node) { left = node; }
	void setRight(KdTree* node) { right = node; }
    void updateBound(BoundingBox b) {
        if (bound.isEmpty()) {
            bound = b;
        }
        else {
            auto bMin = b.getMin();
            auto bMax = b.getMax();
            auto boundMin = bound.getMin();
            auto minX = std::min(boundMin.x, bMin.x);
            auto minY = std::min(boundMin.y, bMin.y);
            auto minZ = std::min(boundMin.z, bMin.z);
            bound.setMin(glm::dvec3(minX, minY, minZ));
            auto boundMax = bound.getMax();
            auto maxX = std::max(boundMax.x, bMax.x);
            auto maxY = std::max(boundMax.y, bMax.y);
            auto maxZ = std::max(boundMax.z, bMax.z);
            bound.setMax(glm::dvec3(maxX, maxY, maxZ));
            // printf("new bound (%f, %f, %f) (%f, %f, %f)\n", minX, minY, minZ, maxX, maxY, maxZ);
        }
    }
    void updateBound(glm::dvec3 bMin, glm::dvec3 bMax) {
        if (bound.isEmpty()) {
            bound = BoundingBox(bMin, bMax);
        }
        else {
            auto boundMin = bound.getMin();
            auto minX = std::min(boundMin.x, bMin.x);
            auto minY = std::min(boundMin.y, bMin.y);
            auto minZ = std::min(boundMin.z, bMin.z);
            bound.setMin(glm::dvec3(minX, minY, minZ));
            auto boundMax = bound.getMax();
            auto maxX = std::max(boundMax.x, bMax.x);
            auto maxY = std::max(boundMax.y, bMax.y);
            auto maxZ = std::max(boundMax.z, bMax.z);
            bound.setMax(glm::dvec3(maxX, maxY, maxZ));
        }
    }
	void addObject(Obj* object) { objects.emplace_back(object);}
	auto& getObjects() { return objects; }
    int getNumObjects() { return objects.size(); }
    bool isLeaf() { return !left && !right; }
};