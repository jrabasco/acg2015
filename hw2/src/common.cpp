#include <nori/object.h>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <boost/math/special_functions/fpclassify.hpp>

#if defined(PLATFORM_LINUX)
#include <malloc.h>
#endif

#if defined(PLATFORM_WINDOWS)
#include <windows.h>
#endif

#if defined(PLATFORM_MACOS)
#include <sys/sysctl.h>
#endif

#if !defined(L1_CACHE_LINE_SIZE)
#define L1_CACHE_LINE_SIZE 64
#endif

NORI_NAMESPACE_BEGIN

Color3f Color3f::toSRGB() const {
	Color3f result;

	for (int i=0; i<3; ++i) {
		float value = coeff(i);

		if (value <= 0.0031308f)
			result[i] = 12.92f * value;
		else
			result[i] = (1.0f + 0.055f) 
				* std::pow(value, 1.0f/2.4f) -  0.055f;
	}

	return result;
}

Color3f Color3f::toLinearRGB() const {
	Color3f result;

	for (int i=0; i<3; ++i) {
		float value = coeff(i);

		if (value <= 0.04045f)
			result[i] = value * (1.0f / 12.92f);
		else
			result[i] = std::pow((value + 0.055f)
				* (1.0f / 1.055f), 2.4f);
	}

	return result;
}

bool Color3f::isValid() const {
	for (int i=0; i<3; ++i) {
		float value = coeff(i);
		int cl = boost::math::fpclassify(value);
		if (value < 0 || cl == FP_INFINITE || cl == FP_NAN)
			return false;
	}
	return true;
}


float Color3f::getLuminance() const {
	return coeff(0) * 0.212671f + coeff(1) * 0.715160f + coeff(2) * 0.072169f;
}

Transform::Transform(const Eigen::Matrix4f &trafo) 
	: m_transform(trafo), m_inverse(trafo.inverse()) { }

QString Transform::toString() const {
	std::ostringstream oss;
	oss << m_transform.format(Eigen::IOFormat(4, 0, ", ", ";\n", "", "", "[", "]"));
	return QString(oss.str().c_str());
}
QString Transform::toLineString() const {
	std::ostringstream oss;
        for(int row = 0; row < 4; ++row){
            if(row > 0) oss << "; ";
            for(int col = 0; col < 4; ++col){
                if(col > 0) oss << ", ";
                oss << m_transform(row, col);
            }
        }
	return QString(oss.str().c_str());
}

Vector3f squareToUniformSphere(const Point2f &sample) {
	float z = 1.0f - 2.0f * sample.y();
	float r = std::sqrt(std::max((float) 0.0f, 1.0f - z*z));
	float sinPhi, cosPhi;
	sincosf(2.0f * M_PI * sample.x(), &sinPhi, &cosPhi);
	return Vector3f(r * cosPhi, r * sinPhi, z);
}

Vector3f squareToUniformHemisphere(const Point2f &sample) {
	float cosTheta = sample.x();
	float sinTheta = std::sqrt(std::max((float) 0, 1-cosTheta*cosTheta));

	float sinPhi, cosPhi;
	sincosf(2.0f * M_PI * sample.y(), &sinPhi, &cosPhi);

	return Vector3f(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
}

Point2f squareToUniformDisk(const Point2f &sample) {
	float r = std::sqrt(sample.x());
	float sinPhi, cosPhi;
	sincosf(2.0f * M_PI * sample.y(), &sinPhi, &cosPhi);

	return Point2f(
		cosPhi * r,
		sinPhi * r
	);
}

Point2f squareToUniformDiskConcentric(const Point2f &sample) {
	float r1 = 2.0f*sample.x() - 1.0f;
	float r2 = 2.0f*sample.y() - 1.0f;

	Point2f coords;
	if (r1 == 0 && r2 == 0) {
		coords = Point2f(0, 0);
	} else if (r1 > -r2) { /* Regions 1/2 */
		if (r1 > r2)
			coords = Point2f(r1, (M_PI/4.0f) * r2/r1);
		else
			coords = Point2f(r2, (M_PI/4.0f) * (2.0f - r1/r2));
	} else { /* Regions 3/4 */
		if (r1<r2)
			coords = Point2f(-r1, (M_PI/4.0f) * (4.0f + r2/r1));
		else 
			coords = Point2f(-r2, (M_PI/4.0f) * (6.0f - r1/r2));
	}

	Point2f result;
	sincosf(coords.y(), &result[1], &result[0]);
	return result*coords.x();
}

Point2f squareToUniformTriangle(const Point2f &sample) {
	float a = std::sqrt(1.0f - sample.x());
	return Point2f(1 - a, a * sample.y());
}

float intervalToTent(float sample) {
	float sign;

	if (sample < 0.5f) {
		sign = 1;
		sample *= 2;
	} else {
		sign = -1;
		sample = 2 * (sample - 0.5f);
	}

	return sign * (1 - std::sqrt(sample));
}

Point2f squareToTent(const Point2f &sample) {
	return Point2f(
		intervalToTent(sample.x()),
		intervalToTent(sample.y())
	);
}

Point2f squareToTriangle(const Point2f &sample) {
	float a = std::sqrt(1.0f - sample.x());
	return Point2f(1 - a, a * sample.y());
}

Vector3f sphericalDirection(float theta, float phi) {
	float sinTheta, cosTheta, sinPhi, cosPhi;

	sincosf(theta, &sinTheta, &cosTheta);
	sincosf(phi, &sinPhi, &cosPhi);

	return Vector3f(
		sinTheta * cosPhi,
		sinTheta * sinPhi,
		cosTheta
	);
}

Point2f sphericalCoordinates(const Vector3f &v) {
	Point2f result(
		std::acos(v.z()),
		std::atan2(v.y(), v.x())
	);
	if (result.y() < 0)
		result.y() += 2*M_PI;
	return result;
}

void coordinateSystem(const Vector3f &a, Vector3f &b, Vector3f &c) {
	if (std::abs(a.x()) > std::abs(a.y())) {
		float invLen = 1.0f / std::sqrt(a.x() * a.x() + a.z() * a.z());
		c = Vector3f(a.z() * invLen, 0.0f, -a.x() * invLen);
	} else {
		float invLen = 1.0f / std::sqrt(a.y() * a.y() + a.z() * a.z());
		c = Vector3f(0.0f, a.z() * invLen, -a.y() * invLen);
	}
	b = c.cross(a);
}

void *allocAligned(size_t size) {
#if defined(PLATFORM_WINDOWS)
	return _aligned_malloc(size, L1_CACHE_LINE_SIZE);
#elif defined(PLATFORM_MACOS)
	/* OSX malloc already returns 16-byte aligned data suitable
	   for AltiVec and SSE computations */
	return malloc(size);
#else
	return memalign(L1_CACHE_LINE_SIZE, size);
#endif
}

void freeAligned(void *ptr) {
#if defined(PLATFORM_WINDOWS)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

int getCoreCount() {
#if defined(PLATFORM_WINDOWS)
	SYSTEM_INFO sys_info;
	GetSystemInfo(&sys_info);
	return sys_info.dwNumberOfProcessors;
#elif defined(PLATFORM_MACOS)
	int nprocs;
	size_t nprocsSize = sizeof(int);
	if (sysctlbyname("hw.activecpu", &nprocs, &nprocsSize, NULL, 0))
		throw NoriException("Could not detect the number of processors!");
	return (int) nprocs;
#else
	return sysconf(_SC_NPROCESSORS_CONF);
#endif
}

QString indent(const QString &string, int amount) {
	QString result = string;
	result.replace("\n", QString("\n") + QString(" ").repeated(amount));
	return result;
}

float fresnel(float cosThetaI, float extIOR, float intIOR) {
	float etaI = extIOR, etaT = intIOR;

	if (extIOR == intIOR)
		return 0.0f;

	/* Swap the indices of refraction if the interaction starts
	   at the inside of the object */
	if (cosThetaI < 0.0f) {
		std::swap(etaI, etaT);
		cosThetaI = -cosThetaI;
	}

	/* Using Snell's law, calculate the squared sine of the
	   angle between the normal and the transmitted ray */
	float eta = etaI / etaT,
		  sinThetaTSqr = eta*eta * (1-cosThetaI*cosThetaI);

	if (sinThetaTSqr > 1.0f)
		return 1.0f;  /* Total internal reflection! */

	float cosThetaT = std::sqrt(1.0f - sinThetaTSqr);

	float Rs = (etaI * cosThetaI - etaT * cosThetaT)
	         / (etaI * cosThetaI + etaT * cosThetaT);
	float Rp = (etaT * cosThetaI - etaI * cosThetaT)
	         / (etaT * cosThetaI + etaI * cosThetaT);

	return (Rs * Rs + Rp * Rp) / 2.0f;
}

NORI_NAMESPACE_END
