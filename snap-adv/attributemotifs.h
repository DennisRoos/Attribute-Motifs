#ifndef snap_attributemotifs_h
#define snap_attributemotifs_h

#include "Snap.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <regex>

typedef TKeyDat<TInt, TInt> TIntPair;
typedef TKeyDat<TInt, TIntPr> TIntTriple;
typedef TKeyDat<TInt, TIntIntPrPr> TIntQuadruple;

#ifndef snap_temporalmotifs_h



// Simple one-dimensional, two-dimensional, and three-dimensional counter
// classes with default intialization.
class Counter1D {
public:
	Counter1D(int m = 0) : m_(m) {
		if (m > 0) {
			data_ = TUInt64V(m);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i) const { return data_[i]; }
	TUInt64& operator()(int i) { return data_[i]; }
	int m() { return m_; }

private:
	int m_;
	TUInt64V data_;
};

class Counter2D {
public:
	Counter2D(int m = 0, int n = 0) : m_(m), n_(n) {
		if (m * n > 0) {
			data_ = TUInt64V(m * n);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i, int j) const { return data_[i + j * m_]; }
	TUInt64& operator()(int i, int j) { return data_[i + j * m_]; }
	int m() { return m_; }
	int n() { return n_; }

private:
	int m_;
	int n_;
	TUInt64V data_;
};

class Counter3D {
public:
	Counter3D(int m = 0, int n = 0, int p = 0) : m_(m), n_(n), p_(p) {
		if (m * n * p > 0) {
			data_ = TUInt64V(m * n * p);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i, int j, int k) const {
		return data_[i + j * m_ + k * m_ * n_];
	}
	TUInt64& operator()(int i, int j, int k) {
		return data_[i + j * m_ + k * m_ * n_];
	}
	int m() { return m_; }
	int n() { return n_; }
	int p() { return p_; }

private:
	int m_;
	int n_;
	int p_;
	TUInt64V data_;
};

#endif

class Counter4D {
public:
	Counter4D(int m = 0, int n = 0, int p = 0, int q = 0) : m_(m), n_(n), p_(p), q_(q) {
		if (m * n * p * q > 0) {
			data_ = TUInt64V(m * n * p * q);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i, int j, int k, int l) const {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_];
	}
	TUInt64& operator()(int i, int j, int k, int l) {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_];
	}
	int m() { return m_; }
	int n() { return n_; }
	int p() { return p_; }
	int q() { return q_; }

private:
	int m_;
	int n_;
	int p_;
	int q_;
	TUInt64V data_;
};

class Counter5D {
public:
	Counter5D(int m = 0, int n = 0, int p = 0, int q = 0, int r = 0) : m_(m), n_(n), p_(p), q_(q), r_(r)
	{
		if (m * n * p * q * r > 0) {
			data_ = TUInt64V(m * n * p * q * r);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i, int j, int k, int l, int m) const {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_];
	}
	TUInt64& operator()(int i, int j, int k, int l, int m) {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_];
	}
	int m() { return m_; }
	int n() { return n_; }
	int p() { return p_; }
	int q() { return q_; }
	int r() { return r_; }

private:
	int m_;
	int n_;
	int p_;
	int q_;
	int r_;
	TUInt64V data_;
};

class Counter6D {
public:
	Counter6D(int m = 0, int n = 0, int p = 0, int q = 0, int r = 0, int s = 0) : m_(m), n_(n), p_(p), q_(q), r_(r), s_(s)
	{
		if (m * n * p * q * r * s > 0) {
			data_ = TUInt64V(m * n * p * q * r * s);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i, int j, int k, int l, int m, int n) const {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_ + n * m_ * n_ * p_ * q_ * r_];
	}
	TUInt64& operator()(int i, int j, int k, int l, int m, int n) {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_ + n * m_ * n_ * p_ * q_ * r_];
	}
	int m() { return m_; }
	int n() { return n_; }
	int p() { return p_; }
	int q() { return q_; }
	int r() { return r_; }
	int s() { return s_; }

private:
	int m_;
	int n_;
	int p_;
	int q_;
	int r_;
	int s_;
	TUInt64V data_;
};


class Counter7D {
public:
	Counter7D(int m = 0, int n = 0, int p = 0, int q = 0, int r = 0, int s = 0, int t = 0) : m_(m), n_(n), p_(p), q_(q), r_(r), s_(s), t_(t)
	{
		if (m * n * p * q * r * s * t > 0) {
			data_ = TUInt64V(m * n * p * q * r * s);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i, int j, int k, int l, int m, int n, int o) const {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_ + n * m_ * n_ * p_ * q_ * r_ + o * m_ * n_ * p_ * q_ * r_ * s_];
	}
	TUInt64& operator()(int i, int j, int k, int l, int m, int n, int o) {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_ + n * m_ * n_ * p_ * q_ * r_ + o * m_ * n_ * p_ * q_ * r_ * s_];
	}
	int m() { return m_; }
	int n() { return n_; }
	int p() { return p_; }
	int q() { return q_; }
	int r() { return r_; }
	int s() { return s_; }
	int t() { return t_; }

private:
	int m_;
	int n_;
	int p_;
	int q_;
	int r_;
	int s_;
	int t_;
	TUInt64V data_;
};

class Counter8D {
public:
	Counter8D(int m = 0, int n = 0, int p = 0, int q = 0, int r = 0, int s = 0, int t = 0, int u = 0) : m_(m), n_(n), p_(p), q_(q), r_(r), s_(s), t_(t), u_(u)
	{
		if (m * n * p * q * r * s * t * u > 0) {
			data_ = TUInt64V(m * n * p * q * r * s * t * u);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i, int j, int k, int l, int m, int n, int o, int p) const {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_ + n * m_ * n_ * p_ * q_ * r_ + o * m_ * n_ * p_ * q_ * r_ * s_ + p * m_ * n_ * p_ * q_ * r_ * s_ * t_];
	}
	TUInt64& operator()(int i, int j, int k, int l, int m, int n, int o, int p) {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_ + n * m_ * n_ * p_ * q_ * r_ + o * m_ * n_ * p_ * q_ * r_ * s_ + p * m_ * n_ * p_ * q_ * r_ * s_ * t_];
	}
	int m() { return m_; }
	int n() { return n_; }
	int p() { return p_; }
	int q() { return q_; }
	int r() { return r_; }
	int s() { return s_; }
	int t() { return t_; }
	int u() { return u_; }

private:
	int m_;
	int n_;
	int p_;
	int q_;
	int r_;
	int s_;
	int t_;
	int u_;
	TUInt64V data_;
};


class Counter9D {
public:
	Counter9D(int m = 0, int n = 0, int p = 0, int q = 0, int r = 0, int s = 0, int t = 0, int u = 0, int v = 0) : m_(m), n_(n), p_(p), q_(q), r_(r), s_(s), t_(t), u_(u), v_(v)
	{
		if (m * n * p * q * r * s * t * u * v > 0) {
			data_ = TUInt64V(m * n * p * q * r * s * t * u * v);
			data_.PutAll(0);
		}
	}
	const TUInt64& operator()(int i, int j, int k, int l, int m, int n, int o, int p, int q) const {//fix this one
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_ + n * m_ * n_ * p_ * q_ * r_ + o * m_ * n_ * p_ * q_ * r_ * s_ + p * m_ * n_ * p_ * q_ * r_ * s_ * t_ + q * m_ * n_ * p_ * q_ * r_ * s_ * t_ * u_];
	}
	TUInt64& operator()(int i, int j, int k, int l, int m, int n, int o, int p, int q) {
		return data_[i + j * m_ + k * m_ * n_ + l * m_ * n_ * p_ + m * m_ * n_ * p_ * q_ + n * m_ * n_ * p_ * q_ * r_ + o * m_ * n_ * p_ * q_ * r_ * s_ + p * m_ * n_ * p_ * q_ * r_ * s_ * t_ + q * m_ * n_ * p_ * q_ * r_ * s_ * t_ * u_];
	}
	int m() { return m_; }
	int n() { return n_; }
	int p() { return p_; }
	int q() { return q_; }
	int r() { return r_; }
	int s() { return s_; }
	int t() { return t_; }
	int u() { return u_; }
	int v() { return v_; }

private:
	int m_;
	int n_;
	int p_;
	int q_;
	int r_;
	int s_;
	int t_;
	int u_;
	int v_;
	TUInt64V data_;
};


// Triad edge data consists of a neighbor, a direction, an indicator of whether
// the edge connects with wich endpoint (u or v) and a layer.
class MultTriadEdgeData {
public:
	MultTriadEdgeData() {}
	MultTriadEdgeData(int _nbr, int _dir, int _u_or_v, int _l, int a) :
		nbr(_nbr), dir(_dir), u_or_v(_u_or_v), l(_l), att(a) {}
	int nbr;     // Which neighbor of the center node
	int dir;     // Outgoing (0) or incoming (1) direction
	int u_or_v;  // Points to first end point u (0) or second end point v (1)
	int l;       // Layer
	int att;	 // attribute of the neighbor
};

// Star edge data consists of a neighbor, a direction and a layer.
class MultStarEdgeData {
public:
	MultStarEdgeData() {}
	MultStarEdgeData(int _nbr, int _dir, int _l, int _a) : nbr(_nbr), dir(_dir), l(_l), att(_a) {}
	int nbr;  // Which neighbor of the center node
	int dir;  // Outgoing (0) or incoming (1) direction
	int l;    // Layer
	int att;  // Attribute of the neighbor node
};

// Main temporal motif counting class.  This implementation has support for
// counting motifs with three multilayer temporal edges on two or three nodes.  
class MultTempMotifCounter {
public:
	// Reads directed multilayer temporal graph data from the specified file, which must have
	// the following format:
	//    source_node destination_node unix_timestamp layer_id
	// Note that both timestamps and layer_id's must be positive integers
	// preferably the layer_id's should be 0,1,...,nrlayers-1 without gaps,
	// otherwise non-existing layers will be assumed to exist.
	// Timestamps of -1 are allowed to indicate no known timestamp!
	// Note that no correct order is guaranteed for such edges.
	// If sim is set to true, all layers will be considered one single layer
	MultTempMotifCounter(const TStr& filename, bool sim, bool att);

	// Count all three multilayer temporal edge, two-node delta-temporal lambda-multilayer motifs and fills the
	// counter counts with the results.
	void Count3MTEdge2Node(double delta, Counter8D& counts);

	// Similar to Count3MTEdge2Node() except only counts motif instances
	// for a given pair of nodes u and v and specifies the source and destination
	// node.  The counts format (disregarding layers) is:
	//   counts(0, 0, 0): u --> v, u --> v, u --> v
	//   counts(1, 1, 1): v --> u, v --> u, v --> u
	//   counts(0, 1, 1): u --> v, v --> u, v --> u
	//   counts(1, 0, 0): v --> u, u --> v, u --> v
	//   counts(0, 1, 0): u --> v, v --> u, u --> v
	//   counts(1, 0, 1): v --> u, u --> v, v --> u
	//   counts(0, 0, 1): u --> v, u --> v, v --> u
	//   counts(1, 1, 0): v --> u, v --> u, u --> v
	void Count3MTEdge2Node(int u, int v, double delta, Counter8D& counts);

	// Counts 3-edge, 3-node star motifs and places the results in pre_counts,
	// pos_counts, and mid_counts.  Counts take the following structure (with
	// center node c):
	//
	//     pre: {c, u}, {c, u}, {c, v}
	//     pos: {c, u}, {c, v}, {c, v}
	//     mid: {c, u}, {c, v}, {c, u}
	//
	// The indices in the counter correspond to the direction with the center
	// node as the reference (0 indexes outgoing and 1 indexes incoming edges to
	// the center).  For example,
	//
	//     pos_counts(0, x, 1, x, 0, x): c --> u, v --> c, c --> v, with the x's denoting their layers
	void Count3MTEdge3NodeStars(double delta, Counter9D& pre_counts,
		Counter9D& pos_counts, Counter9D& mid_counts);

	// Counts 3-edge triad events and places the result in counts (disregarding layers):
	//
	//    counts(0, 0, 0): u --> v, w --> v, u --> w (M_{1,3})
	//    counts(0, 0, 1): u --> v, w --> v, w --> u (M_{1,4})
	//    counts(0, 1, 0): u --> v, v --> w, u --> w (M_{2,3})
	//    counts(0, 1, 1): u --> v, v --> w, w --> u (M_{2,4})
	//    counts(1, 0, 0): u --> v, w --> u, v --> w (M_{3,5})
	//    counts(1, 0, 1): u --> v, w --> u, w --> v (M_{3,6})
	//    counts(1, 1, 0): u --> v, u --> w, v --> w (M_{4,5})
	//    counts(1, 1, 1): u --> v, u --> w, w --> v (M_{4,6})
	void Count3MTEdgeTriads(double delta, Counter9D& counts);

	// Counts all 3-edge, {2,3}-node multilayer temporal motifs and places the result in
	// counts such that counts(i, j, k) corresponds to motif M_{i,j,k}.
	void Count3MTEdge23Node(double delta, Counter3D& counts2, Counter4D& counts3);


	int nodecount;

	TInt nrlayers_;

	TInt nratts_;
	std::vector<int> attribute_data;
private:

	// Get all triangles in the static graph, (Us(i), Vs(i), Ws(i)) is the ith
	// triangle.
	void GetAllStaticTriangles(TIntV& Us, TIntV& Vs, TIntV& Ws);
	// Fills nbrs with all neighbours (ignoring direction) of the node in the
	// static graph.
	void GetAllNeighbors(int node, TIntV& nbrs);
	// Fills nodes with a vector of all nodes in the static graph.
	void GetAllNodes(TIntV& nodes);

	// Checks whether or not there is an edge along the static edge (u, v)
	bool HasEdges(int u, int v);

	// A simple wrapper for adding triad edge data
	void AddMultTriadEdgeData(TVec<MultTriadEdgeData>& events, TVec<TIntPair>& ts_indices,
		int& index, int u, int v, int nbr, int key1, int key2, int att);
	// A simple wrapper for adding star edge data  
	void AddMultStarEdgeData(TVec<TIntPair>& ts_indices, TVec<MultStarEdgeData>& events,
		int& index, int u, int v, int nbr, int key, int att);
	// Another simple wrapper for adding star edge data
	void AddStarEdges(TVec<TIntQuadruple>& combined, int& index, int u, int v, int key);

	// Ensure that edges with equal timestamps remain in the same order as they are initially added to the sequence.
	void FixSameTimestampOrder(TVec<TIntPair>& order);
	void FixSameTimestampOrder(TVec<TIntQuadruple>& order);

	// Directed graph from ignoring timestamps
	PNGraph static_graph_;

	// Core data structure for storing temporal edges.  temporal_data_[u](v) is a
	// list of temporal edges along the static edge (u, v).
	TVec< THash<TInt, TIntV> > temporal_data_;

	// Core data structure for storing 
	TVec< THash<TInt, TIntV> > layer_data_;

	// Indicates the highest layer label, assuming no gaps it also indicates the number of layers

	// Indicates whether partial timing was detected and chooses most efficient counting scheme accordingly
	TBool partial_;
};

// This class exhaustively counts all size^3 three-edge multilayer temporal motifs in an
// alphabet of a given size.
class ThreeMTEdgeMotifCounter {
public:

	// Initialize counter with a given alphabet size and number of layers
	ThreeMTEdgeMotifCounter(int size, int nrlayers, int nratts) : size_(size), nrlayers_(nrlayers), nratts_(nratts) {}

	// Count all three-edge motifs with corresponding timestamps.  Each integer in
	// the event_string must belong to the set {0, 1, ..., size - 1}.  The function
	// stores the results in the counter, where counts(e, f, g) is the motif consisting
	// of the ordered edges e, f, g.
	void Count(const TIntV& event_string, const TIntV& timestamps, const TIntV& layers, const TInt& attrSrc, const TInt& attrDst,
		double delta, Counter8D& counts);

private:
	void IncrementCounts(int event, int layer, const TInt& attSrc, const TInt& attDst);
	void DecrementCounts(int event, int layer);
	Counter2D counts1_;
	Counter4D counts2_;
	Counter8D counts3_;
	int size_;  // alphabet size
	int nrlayers_; // number of layers
	int nratts_;
};

//***********************************************************************************************

// This class exhaustively counts all size^3 three-edge multilayer temporal motifs in an
// alphabet of a given size.
class ThreeMpTEdgeMotifCounter {
public:
	// Initialize counter with a given alphabet size and number of layers
	ThreeMpTEdgeMotifCounter(int size, int nrlayers, int nratts) : size_(size), nrlayers_(nrlayers), nratts_(nratts) {}

	// Count all three-edge motifs with corresponding timestamps.  Each integer in
	// the event_string must belong to the set {0, 1, ..., size - 1}.  The function
	// stores the results in the counter, where counts(e, f, g) is the motif consisting
	// of the ordered edges e, f, g.
	void Count(const TIntV& event_string, const TIntV& timestamps, const TIntV& layers, const TInt& attrSrc, const TInt& attrDst,
		double delta, Counter8D& counts);

private:
	void PreProcess(int event, int layer, int att1, int att2);
	void IncrementCounts(int event, int layer, int att1, int att2);
	void DecrementCounts(int event, int layer, int att1, int att2);
	Counter2D counts1_;
	Counter4D counts2_;
	Counter8D counts3_;
	Counter2D pcounts1_;
	Counter4D pcounts2_;
	int size_;  // alphabet size
	int nrlayers_; // number of layers
	int nratts_;//number of attributes
};

//***********************************************************************************************

// Base class for 3-edge, 3-node star and triangle counters.  The template type
// describes the data needed when processing an edge.
template <typename EdgeData>
class StarTriad3MTEdgeCounter {
public:
	StarTriad3MTEdgeCounter() {}
	void Count(const TVec<EdgeData>& events, const TIntV& timestamps, double delta);

protected:
	// These methods depend on the motif type (star or triad).
	virtual void InitializeCounters() = 0;
	virtual void PopPre(const EdgeData& event) = 0;
	virtual void PopPos(const EdgeData& event) = 0;
	virtual void PushPre(const EdgeData& event) = 0;
	virtual void PushPos(const EdgeData& event) = 0;
	virtual void ProcessCurrent(const EdgeData& event) = 0;
};

//***********************************************************************************************

// Base class for 3-edge, 3-node star and triangle counters.  The template type
// describes the data needed when processing an edge. Used for partially timed edge sequences.
template <typename EdgeData>
class StarTriad3MpTEdgeCounter {
public:
	StarTriad3MpTEdgeCounter() {}
	void Count(const TVec<EdgeData>& events, const TIntV& timestamps, double delta);

protected:
	// These methods depend on the motif type (star or triad).
	virtual void InitializeCounters() = 0;
	virtual void PopPre(const EdgeData& event) = 0;
	virtual void PopPos(const EdgeData& event) = 0;
	virtual void PushPre(const EdgeData& event) = 0;
	virtual void PushPos(const EdgeData& event) = 0;
	virtual void ProcessCurrent(const EdgeData& event) = 0;

	virtual void PreProcPopPos(const EdgeData& event) = 0;
	virtual void PreProcPushPre(const EdgeData& event) = 0;
	virtual void PreProcPushPos(const EdgeData& event) = 0;
	virtual void PreProcProcessCurrent(const EdgeData& event) = 0;
};

//***********************************************************************************************

// Class for counting star motifs with a given center node.
class ThreeMTEdgeStarCounter : public StarTriad3MTEdgeCounter<MultStarEdgeData> {
public:
	// Construct class with maximum number of neighbor nodes.  Each processed edge
	// consists of a neighbour and a direction where the neighbour is represented by
	// an int belong to the set {0, 1, ..., max_nodes - 1}.
	ThreeMTEdgeStarCounter(int max_nodes, int nrlayers, int nratts, int center_att) : max_nodes_(max_nodes), nrlayers_(nrlayers), nratts_(nratts), c_att(center_att) {}

	// Counting conventions follow MultTempMotifCounter::Count3MTEdge3NodeStars().
	int PreCount(int dir1, int lay1, int dir2, int lay2, int dir3, int lay3, int att1, int att2, int att3) {
		return pre_counts_(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
	}
	int PosCount(int dir1, int lay1, int dir2, int lay2, int dir3, int lay3, int att1, int att2, int att3) {
		return pos_counts_(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
	}
	int MidCount(int dir1, int lay1, int dir2, int lay2, int dir3, int lay3, int att1, int att2, int att3) {
		return mid_counts_(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
	}

protected:
	void InitializeCounters();
	void PopPre(const MultStarEdgeData& event);
	void PopPos(const MultStarEdgeData& event);
	void PushPre(const MultStarEdgeData& event);
	void PushPos(const MultStarEdgeData& event);
	void ProcessCurrent(const MultStarEdgeData& event);

private:
	int max_nodes_;
	int nrlayers_;
	int nratts_;
	int c_att; 
	Counter5D pre_sum_;
	Counter5D pos_sum_;
	Counter5D mid_sum_;
	Counter9D pre_counts_;
	Counter9D pos_counts_;
	Counter9D mid_counts_;
	Counter4D pre_nodes_;
	Counter4D pos_nodes_;
};

//***********************************************************************************************

// Class for counting star motifs with a given center node, given a partially timed sequence of edges.
class ThreeMpTEdgeStarCounter : public StarTriad3MpTEdgeCounter<MultStarEdgeData> {
public:
	// Construct class with maximum number of neighbor nodes.  Each processed edge
	// consists of a neighbor and a direction where the neighbor is represented by
	// an int belong to the set {0, 1, ..., max_nodes - 1}.
	ThreeMpTEdgeStarCounter(int max_nodes, int nrlayers, int nratts, int center_att) : max_nodes_(max_nodes), nrlayers_(nrlayers), nratts_(nratts), c_att(center_att) {}

	// Counting conventions follow MultTempMotifCounter::Count3MTEdge3NodeStars().
	int PreCount(int dir1, int lay1, int dir2, int lay2, int dir3, int lay3, int att1, int att2, int att3) {
		return pre_counts_(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
	}
	int PosCount(int dir1, int lay1, int dir2, int lay2, int dir3, int lay3, int att1, int att2, int att3) {
		return pos_counts_(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
	}
	int MidCount(int dir1, int lay1, int dir2, int lay2, int dir3, int lay3, int att1, int att2, int att3) {
		return mid_counts_(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
	}

protected:
	void InitializeCounters();
	void PopPre(const MultStarEdgeData& event);
	void PopPos(const MultStarEdgeData& event);
	void PushPre(const MultStarEdgeData& event);
	void PushPos(const MultStarEdgeData& event);
	void ProcessCurrent(const MultStarEdgeData& event);

	void PreProcPopPos(const MultStarEdgeData& event);
	void PreProcPushPre(const MultStarEdgeData& event);
	void PreProcPushPos(const MultStarEdgeData& event);
	void PreProcProcessCurrent(const MultStarEdgeData& event);

private:
	int max_nodes_;
	int nrlayers_;
	int nratts_;
	int c_att;

	Counter5D pre_sum_;
	Counter5D pos_sum_;
	Counter5D mid_sum_;
	Counter9D pre_counts_;
	Counter9D pos_counts_;
	Counter9D mid_counts_;
	Counter3D pre_nodes_;
	Counter3D pos_nodes_;

	Counter4D p_pre_mid_;

	Counter4D p_pre_sum_;
	Counter4D p_pos_sum_;
	Counter4D p_mid_sum_;
	Counter3D p_pre_nodes_;
	Counter3D p_pos_nodes_;
};

//***********************************************************************************************

// Class for counting triangle motifs that contain a specific undirected edge.
class ThreeMTEdgeTriadCounter : public StarTriad3MTEdgeCounter<MultTriadEdgeData> {
public:
	// Construct class with maximum number of neighbor nodes.  Each processed edge
	// consists of a neighbor, a direction, and an indicator of which end point it
	// connects with.  Each neighbor is represented by an int belong to the set
	// {0, 1, ..., max_nodes - 1}.
	ThreeMTEdgeTriadCounter(int max_nodes, int node_u, int node_v, int nrlayers, int nratts, std::vector<int> attributes, int uatt, int vatt) :
		max_nodes_(max_nodes), node_u_(node_u), node_v_(node_v), nrlayers_(nrlayers), nratts_(nratts), attribute_data(attributes), u_att(uatt), v_att(vatt) {}

	// Counting conventions follow MultTempMotifCounter::Count3MTEdgeTriads().
	int Counts(int dir1, int lay1, int dir2, int lay2, int dir3, int lay3, int at1, int at2, int at3) {
		return triad_counts_(dir1, lay1, dir2, lay2, dir3, lay3, at1, at2, at3);
	}
	std::vector<int> attribute_data;


protected:
	void InitializeCounters();
	void PopPre(const MultTriadEdgeData& event);
	void PopPos(const MultTriadEdgeData& event);
	void PushPre(const MultTriadEdgeData& event);
	void PushPos(const MultTriadEdgeData& event);
	void ProcessCurrent(const MultTriadEdgeData& event);
	bool IsEdgeNode(int nbr) { return nbr == node_u_ || nbr == node_v_; }

private:
	int max_nodes_;
	Counter6D pre_sum_;
	Counter6D pos_sum_;
	Counter6D mid_sum_;
	Counter9D triad_counts_;
	Counter4D pre_nodes_;
	Counter4D pos_nodes_;
	// Two end points of the edge whose triangles this class counts.
	int node_u_;
	int node_v_;
	int u_att;
	int v_att;
	int nrlayers_;
	int nratts_;
};


//***********************************************************************************************

// Class for counting triangle motifs that contain a specific undirected edge, given a partially timed sequence of edges.
class ThreeMpTEdgeTriadCounter : public StarTriad3MpTEdgeCounter<MultTriadEdgeData> {
public:
	// Construct class with maximum number of neighbor nodes.  Each processed edge
	// consists of a neighbor, a direction, and an indicator of which end point it
	// connects with.  Each neighbor is represented by an int belong to the set
	// {0, 1, ..., max_nodes - 1}.
	ThreeMpTEdgeTriadCounter(int max_nodes, int node_u, int node_v, int nrlayers, int nratts, std::vector<int> attributes, int uatt, int vatt) :
		max_nodes_(max_nodes), node_u_(node_u), node_v_(node_v), nrlayers_(nrlayers), nratts_(nratts), attribute_data(attributes), u_att(uatt), v_att(vatt) {}

	// Counting conventions follow MultTempMotifCounter::Count3MTEdgeTriads().
	int Counts(int dir1, int lay1, int dir2, int lay2, int dir3, int lay3, int at1, int at2, int at3) {
		return triad_counts_(dir1, lay1, dir2, lay2, dir3, lay3, at1, at2, at3);
	}

	std::vector<int> attribute_data;

protected:
	void InitializeCounters();
	void PopPre(const MultTriadEdgeData& event);
	void PopPos(const MultTriadEdgeData& event);
	void PushPre(const MultTriadEdgeData& event);
	void PushPos(const MultTriadEdgeData& event);
	void ProcessCurrent(const MultTriadEdgeData& event);
	bool IsEdgeNode(int nbr) { return nbr == node_u_ || nbr == node_v_; }

	void PreProcPopPos(const MultTriadEdgeData& event);
	void PreProcPushPre(const MultTriadEdgeData& event);
	void PreProcPushPos(const MultTriadEdgeData& event);
	void PreProcProcessCurrent(const MultTriadEdgeData& event);

private:
	int max_nodes_;
	Counter6D pre_sum_;
	Counter6D pos_sum_;
	Counter6D mid_sum_;
	Counter9D triad_counts_;
	Counter4D pre_nodes_;
	Counter4D pos_nodes_;

	Counter5D p_pre_mid_;

	Counter5D p_pre_sum_;
	Counter5D p_pos_sum_;
	Counter5D p_mid_sum_;
	Counter4D p_pre_nodes_;
	Counter4D p_pos_nodes_;
	// Two end points of the edge whose triangles this class counts.
	int node_u_;
	int node_v_;
	int nrlayers_;
	int nratts_;
	int u_att;
	int v_att;
};



#endif  // snap_temporalmotifs_h
