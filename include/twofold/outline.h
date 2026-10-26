#pragma once

#include <cstdint>
#include <vector>

#include "math.h"

namespace twofold {

// 8-bit coverage. Strides are in bytes, so one channel of an interleaved image can be read in place.
struct AlphaImage {
    const uint8_t* pixels = nullptr;
    int width = 0;
    int height = 0;
    int rowStride = 0;
    int pixelStride = 1;
};

// Closed boundaries where alpha crosses threshold, interpolated between pixel centers. Units are
// pixels with the origin at the top-left corner and y down. Pixels outside the image are transparent.
std::vector<std::vector<Vec2>> traceAlpha(const AlphaImage& image, uint8_t threshold = 128);

// Drops vertices while keeping every point of the closed contour within tolerance of the result, which
// never crosses itself. Corners move to where the straight runs beside them meet, which restores
// corners cut by tracing.
std::vector<Vec2> simplifyContour(const std::vector<Vec2>& contour, float tolerance);

// Convex polygons with positive area that tile the region inside the contours. Contours are filled
// even-odd, so holes and the islands inside them may wind either way.
std::vector<std::vector<Vec2>> convexDecompose(const std::vector<std::vector<Vec2>>& contours, int maxVertices = 8);

// traceAlpha, simplifyContour and convexDecompose in one call. The simplified contours never touch.
std::vector<std::vector<Vec2>> alphaPieces(const AlphaImage& image, uint8_t threshold = 128, float tolerance = 0.5f,
                                           int maxVertices = 8);

}
