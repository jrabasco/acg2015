/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2012 by Wenzel Jakob and Steve Marschner.

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/* =======================================================================
     This file contains classes for parallel rendering of "image blocks".
 * ======================================================================= */

#if !defined(__PARALLEL_H)
#define __PARALLEL_H

#include <nori/color.h>
#include <nori/vector.h>
#include <QMutex>
#include <QThread>
#include <QElapsedTimer>

#define NORI_BLOCK_SIZE 32 /* Block size used for parallelization */

NORI_NAMESPACE_BEGIN

/**
 * \brief Weighted pixel storage for a rectangular subregion of an image
 *
 * This class implements storage for a rectangular subregion of a
 * larger image that is being rendered. For each pixel, it records color
 * values along with a weight that specifies the accumulated influence of
 * nearby samples on the pixel (according to the used reconstruction filter).
 *
 * When rendering with filters, the samples in a rectangular
 * region will generally also contribute to pixels just outside of 
 * this region. For that reason, this class also stores information about
 * a small border region around the rectangle, whose size depends on the
 * properties of the reconstruction filter.
 */
class ImageBlock : public Eigen::Array<Color4f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> {
public:
	/**
	 * Create a new image block of the specified maximum size
	 * \param size
	 *     Desired maximum size of the block
	 * \param filter
	 *     Samples will be convolved with the image reconstruction
	 *     filter provided here.
	 */
	ImageBlock(const Vector2i &size, const ReconstructionFilter *filter);
        
        /**
         * Load an image block from a bitmap
         */
        ImageBlock(const Bitmap *image);
	
	/// Release all memory
	~ImageBlock();
	
	/// Configure the offset of the block within the main image
	void setOffset(const Point2i &offset) { m_offset = offset; }

	/// Return the offset of the block within the main image
	inline const Point2i &getOffset() const { return m_offset; }
	
	/// Configure the size of the block within the main image
	void setSize(const Point2i &size) { m_size = size; }

	/// Return the size of the block within the main image
	inline const Vector2i &getSize() const { return m_size; }

	/// Return the border size in pixels
	inline int getBorderSize() const { return m_borderSize; }

	/**
	 * \brief Turn the block into a proper bitmap
	 * 
	 * This entails normalizing all pixels and discarding
	 * the border region.
	 */
	Bitmap *toBitmap() const;

	/// Clear all contents
	void clear() { setConstant(Color4f()); }

	/// Record a sample with the given position and radiance value
	void put(const Point2f &pos, const Color3f &value);

	/**
	 * \brief Merge another image block into this one
	 *
	 * During the merge operation, this function locks 
	 * the destination block using a mutex.
	 */
	void put(ImageBlock &b);

	/// Lock the image block (using an internal mutex)
	inline void lock() const { m_mutex.lock(); }
	
	/// Unlock the image block
	inline void unlock() const { m_mutex.unlock(); }

	/// Return a human-readable string summary
	QString toString() const;
protected:
	Point2i  m_offset;
	Vector2i m_size;
	int m_borderSize;
	float *m_filter, m_filterRadius;
	float *m_weightsX, *m_weightsY;
	float m_lookupFactor;
	mutable QMutex m_mutex;
};

/**
 * \brief Spiraling block generator
 *
 * This class can be used to chop up an image into many small
 * rectangular blocks suitable for parallel rendering. The blocks
 * are ordered in spiraling pattern so that the center is
 * rendered first.
 */
class BlockGenerator {
public:
	/**
	 * \brief Create a block generator with
	 * \param size
	 *      Size of the image that should be split into blocks
	 * \param blockSize
	 *      Maximum size of the individual blocks
	 */
	BlockGenerator(const Vector2i &size, int blockSize);
	
	/**
	 * \brief Return the next block to be rendered
	 *
	 * This function is thread-safe
	 *
	 * \return \c false if there were no more blocks
	 */
	bool next(ImageBlock &block);
protected:
	enum EDirection { ERight = 0, EDown, ELeft, EUp };

	Point2i m_block;
	Vector2i m_numBlocks;
	Vector2i m_size;
	int m_blockSize;
	int m_numSteps;
	int m_blocksLeft;
	int m_stepsLeft;
	int m_direction;
	QMutex m_mutex;
	QElapsedTimer m_timer;
};

/**
 * \brief Render thread
 *
 * This class implements the main rendering logic, which consists of
 * fetching work from a scheduler (in the form of rectangular image
 * blocks to be rendered), processing it, and writing the output
 * to a target buffer.
 */
class BlockRenderThread : public QThread {
public:
	/**
	 * \brief Create a new rendering thread that fetches blocks from
	 * the specified block generator and writes output to a big
	 * \ref ImageBlock instance that represents the entire image
	 */
	BlockRenderThread(const Scene *scene, Sampler *sampler,
		BlockGenerator *blockGenerator, ImageBlock *output);

	/// Release all memory
	virtual ~BlockRenderThread();

	/// Main rendering thread loop
	void run();
private:
	const Scene *m_scene;
	BlockGenerator *m_blockGenerator;
	ImageBlock *m_output;
	Sampler *m_sampler;
};

NORI_NAMESPACE_END

#endif /* __PARALLEL_H */
