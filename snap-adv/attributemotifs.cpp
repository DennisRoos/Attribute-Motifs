#include "Snap.h"
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "attributemotifs.h"

#include <iostream>
#include <fstream>
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <regex>

#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Initialization and helper methods for MultTempMotifCounter
MultTempMotifCounter::MultTempMotifCounter(const TStr& filename, bool sim, bool att) {
	// First load the static graph
	static_graph_ = TSnap::LoadEdgeList<PNGraph>(filename, 0, 1);
	int max_nodes = static_graph_->GetMxNId();
	nodecount = max_nodes;
	temporal_data_ = TVec< THash<TInt, TIntV> >(max_nodes);
	layer_data_ = TVec< THash<TInt, TIntV> >(max_nodes);  

	for (int i = 0; i < max_nodes; i++) {
		attribute_data.push_back(-1);
	}

	int max_attributes = 0;
	int max_layer = 0;
	partial_ = false;
	std::string s = filename();

	// Formulate input File Format:
	//   source_node destination_node timestamp layer source_node_attribute destination_node_attribute
	

	std::ifstream infile;
	infile.open(s);
	int src;
	int dst;
	int tim;
	int lay;
	int at1;
	int at2;
	while (infile >> src) {
		infile >> dst;
		infile >> tim;
		infile >> lay;
		infile >> at1;
		infile >> at2;

		if (sim) { lay = 0; }
		if (att) {
			at1 = 0;
			at2 = 0;
		}
		if (lay > max_layer) { max_layer = lay; }
		if (at1 > max_attributes) { max_attributes = at1; }
		if (at2 > max_attributes) { max_attributes = at2; }

		temporal_data_[src](dst).Add(tim);
		layer_data_[src](dst).Add(lay);
		attribute_data[src] = at1;
		attribute_data[dst] = at2;

	}
	nrlayers_ = max_layer + 1;
	nratts_ = max_attributes + 1;
	// TODO: include check to guarantee no negative timestamps other than -1
}

void MultTempMotifCounter::GetAllNodes(TIntV& nodes) {
	nodes = TIntV();
	for (TNGraph::TNodeI it = static_graph_->BegNI();
		it < static_graph_->EndNI(); it++) {
		nodes.Add(it.GetId());
	}
}

bool MultTempMotifCounter::HasEdges(int u, int v) {
	return temporal_data_[u].IsKey(v);
}

void MultTempMotifCounter::GetAllNeighbors(int node, TIntV& nbrs) {
	nbrs = TIntV();
	TNGraph::TNodeI NI = static_graph_->GetNI(node);
	for (int i = 0; i < NI.GetOutDeg(); i++) { nbrs.Add(NI.GetOutNId(i)); }
	for (int i = 0; i < NI.GetInDeg(); i++) {
		int nbr = NI.GetInNId(i);
		if (!NI.IsOutNId(nbr)) { nbrs.Add(nbr); }
	}
}

void MultTempMotifCounter::GetAllStaticTriangles(TIntV& Us, TIntV& Vs, TIntV& Ws) {
	Us.Clr();
	Vs.Clr();
	Ws.Clr();
	// Get degree ordering of the graph
	int max_nodes = static_graph_->GetMxNId();
	TVec<TIntPair> degrees(max_nodes);
	degrees.PutAll(TIntPair(0, 0));
	// Set the degree of a node to be the number of nodes adjacent to the node in
	// the undirected graph.
	TIntV nodes;
	GetAllNodes(nodes);
#pragma omp parallel for schedule(dynamic)  
	for (int node_id = 0; node_id < nodes.Len(); node_id++) {
		int src = nodes[node_id];
		TIntV nbrs;
		GetAllNeighbors(src, nbrs);
		degrees[src] = TIntPair(nbrs.Len(), src);
	}
	degrees.Sort();
	TIntV order = TIntV(max_nodes);
#pragma omp parallel for schedule(dynamic)  
	for (int i = 0; i < order.Len(); i++) {
		order[degrees[i].Dat] = i;
	}

	// Get triangles centered at a given node where that node is the smallest in
	// the degree ordering.
#pragma omp parallel for schedule(dynamic)  
	for (int node_id = 0; node_id < nodes.Len(); node_id++) {
		int src = nodes[node_id];
		int src_pos = order[src];

		// Get all neighbors who come later in the ordering
		TIntV nbrs;
		GetAllNeighbors(src, nbrs);
		TIntV neighbors_higher;
		for (int i = 0; i < nbrs.Len(); i++) {
			int nbr = nbrs[i];
			if (order[nbr] > src_pos) { neighbors_higher.Add(nbr); }
		}

		for (int ind1 = 0; ind1 < neighbors_higher.Len(); ind1++) {
			for (int ind2 = ind1 + 1; ind2 < neighbors_higher.Len(); ind2++) {
				int dst1 = neighbors_higher[ind1];
				int dst2 = neighbors_higher[ind2];
				// Check for triangle formation
				if (static_graph_->IsEdge(dst1, dst2) || static_graph_->IsEdge(dst2, dst1)) {
#pragma omp critical
					{
						Us.Add(src);
						Vs.Add(dst1);
						Ws.Add(dst2);
					}
				}
			}
		}
	}
}

void MultTempMotifCounter::Count3MTEdge23Node(double delta, Counter3D& counts2,  Counter4D& counts3) {
	// This is simply a wrapper function around the counting methods to produce
	// counts in the same way that they were represented in the paper.  This makes
	// it easy to reproduce results and allow SNAP users to make the same
	// measurements on their temporal network data.
	int total_motifs = (nrlayers_ * nrlayers_ * nrlayers_);	
	int total_attributes = (nratts_ * nratts_ * nratts_);

	int a, b, c;
	int p, q;
	int r;


	Counter8D edge_counts = Counter8D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_);
	Count3MTEdge2Node(delta, edge_counts);


	for (int i = 0; i < total_motifs; i++) {
		c = (i / (nrlayers_ * nrlayers_));
		b = (i / nrlayers_) % nrlayers_;
		a = i % nrlayers_;
		for (int j = 0; j <  nratts_ * nratts_; j++) {
			p = j / nratts_;
			q = j % nratts_;
			counts2(0, i, j) = edge_counts(0, a, 0, b, 0, c, p, q) + edge_counts(1, a, 1, b, 1, c, q, p); // M_{6,1}
			counts2(1, i, j) = edge_counts(1, a, 0, b, 0, c, p, q) + edge_counts(0, a, 1, b, 1, c, q, p); // M_{5,2}
			counts2(2, i, j) = edge_counts(0, a, 1, b, 0, c, p, q) + edge_counts(1, a, 0, b, 1, c, q, p); // M_{5,1}
			counts2(3, i, j) = edge_counts(0, a, 0, b, 1, c, p, q) + edge_counts(1, a, 1, b, 0, c, q, p); // M_{6,2}
		}
	}

	Counter9D pre_counts, pos_counts, mid_counts;	

	Count3MTEdge3NodeStars(delta, pre_counts, pos_counts, mid_counts);
	for (int j = 0; j < total_attributes; j++) {
		p = (j / (nratts_ * nratts_));
		q = (j / nratts_) % nratts_;
		r = j % nratts_;
		for (int i = 0; i < total_motifs; i++) {
			c = (i / (nrlayers_ * nrlayers_));
			b = (i / nrlayers_) % nrlayers_;
			a = i % nrlayers_;
			counts3(3, 5, i, j) = mid_counts(1, a, 1, b, 1, c, p, q, r); // M_{1,1}
			counts3(2, 5, i, j) = mid_counts(1, a, 1, b, 0, c, p, q, r); // M_{1,2}
			counts3(2, 3, i, j) = pos_counts(1, a, 1, b, 0, c, p, q, r); // M_{1,5}
			counts3(3, 3, i, j) = pos_counts(1, a, 1, b, 1, c, p, q, r); // M_{1,6}
			counts3(3, 4, i, j) = mid_counts(1, a, 0, b, 1, c, p, q, r); // M_{2,1}
			counts3(2, 4, i, j) = mid_counts(1, a, 0, b, 0, c, p, q, r); // M_{2,2}
			counts3(0, 3, i, j) = pos_counts(1, a, 0, b, 0, c, p, q, r); // M_{2,5}
			counts3(1, 3, i, j) = pos_counts(1, a, 0, b, 1, c, p, q, r); // M_{2,6}
			counts3(0, 5, i, j) = mid_counts(0, a, 0, b, 1, c, p, q, r); // M_{3,1}
			counts3(1, 5, i, j) = mid_counts(0, a, 1, b, 1, c, p, q, r); // M_{3,2}
			counts3(2, 2, i, j) = pos_counts(0, a, 1, b, 0, c, p, q, r); // M_{3,3}
			counts3(3, 2, i, j) = pos_counts(0, a, 1, b, 1, c, p, q, r); // M_{3,4}
			counts3(0, 4, i, j) = mid_counts(0, a, 0, b, 0, c, p, q, r); // M_{4,1}
			counts3(1, 4, i, j) = mid_counts(0, a, 0, b, 1, c, p, q, r); // M_{4,2}
			counts3(0, 2, i, j) = pos_counts(0, a, 0, b, 0, c, p, q, r); // M_{4,3}
			counts3(1, 2, i, j) = pos_counts(0, a, 0, b, 1, c, p, q, r); // M_{4,4}
			counts3(1, 6, i, j) = pre_counts(0, a, 1, b, 0, c, p, q, r); // M_{5,3}
			counts3(1, 7, i, j) = pre_counts(0, a, 1, b, 1, c, p, q, r); // M_{5,4}
			counts3(2, 6, i, j) = pre_counts(1, a, 0, b, 0, c, p, q, r); // M_{5,5}
			counts3(2, 7, i, j) = pre_counts(1, a, 0, b, 1, c, p, q, r); // M_{5,6}
			counts3(0, 6, i, j) = pre_counts(0, a, 0, b, 0, c, p, q, r); // M_{6,3}
			counts3(0, 7, i, j) = pre_counts(0, a, 0, b, 1, c, p, q, r); // M_{6,4}
			counts3(3, 6, i, j) = pre_counts(1, a, 1, b, 0, c, p, q, r); // M_{6,5}
			counts3(3, 7, i, j) = pre_counts(1, a, 1, b, 1, c, p, q, r); // M_{6,6}
		}
	}

	Counter9D triad_counts;
	Count3MTEdgeTriads(delta, triad_counts);
	for (int j = 0; j < total_attributes; j++) {
		p = (j / (nratts_ * nratts_));
		q = (j / nratts_) % nratts_;
		r = j % nratts_;
		for (int i = 0; i < total_motifs; i++) {
			c = (i / (nrlayers_ * nrlayers_));
			b = (i / nrlayers_) % nrlayers_;
			a = i % nrlayers_;
			counts3(2, 1, i, j) = triad_counts(0, a, 0, b, 0, c, p, q, r); // M_{1,3}
			counts3(3, 1, i, j) = triad_counts(0, a, 0, b, 1, c, p, q, r); // M_{1,4}
			counts3(0, 1, i, j) = triad_counts(0, a, 1, b, 0, c, p, q, r); // M_{2,3}
			counts3(1, 1, i, j) = triad_counts(0, a, 1, b, 1, c, p, q, r); // M_{2,4}
			counts3(1, 0, i, j) = triad_counts(1, a, 0, b, 0, c, p, q, r); // M_{3,5}
			counts3(3, 0, i, j) = triad_counts(1, a, 0, b, 1, c, p, q, r); // M_{3,6}
			counts3(0, 0, i, j) = triad_counts(1, a, 1, b, 0, c, p, q, r); // M_{4,5}
			counts3(2, 0, i, j) = triad_counts(1, a, 1, b, 1, c, p, q, r); // M_{4,6}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Two-node (static edge) counting methods
void MultTempMotifCounter::Count3MTEdge2Node(double delta, Counter8D& counts) {
	// Get a vector of undirected edges (so we can use openmp parallel for over it)
	TVec<TIntPair> undir_edges;


	for (TNGraph::TEdgeI it = static_graph_->BegEI(); it < static_graph_->EndEI(); it++) {
		int src = it.GetSrcNId();
		int dst = it.GetDstNId();
		// Only consider undirected edges
		if (src < dst || (dst < src && !static_graph_->IsEdge(dst, src))) {
			undir_edges.Add(TIntPair(src, dst));
		}
	}

	counts = Counter8D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_);
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < undir_edges.Len(); i++) {
		TIntPair edge = undir_edges[i];
		Counter8D local = Counter8D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_);;



		Count3MTEdge2Node(edge.Key, edge.Dat, delta, local);



#pragma omp critical
		{  // Update with local counts
			for (int dir1 = 0; dir1 < 2; ++dir1) {
				for (int lay1 = 0; lay1 < nrlayers_; ++lay1) {
					for (int dir2 = 0; dir2 < 2; ++dir2) {
						for (int lay2 = 0; lay2 < nrlayers_; ++lay2) {
							for (int dir3 = 0; dir3 < 2; ++dir3) {
								for (int lay3 = 0; lay3 < nrlayers_; ++lay3) {
									for (int att1 = 0; att1 < nratts_; ++att1) {
										for (int att2 = 0; att2 < nratts_; ++att2) {
											counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2) +=
												local(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void MultTempMotifCounter::FixSameTimestampOrder(TVec<TIntPair>& order) {
	int tempDat, i;
	for (int j = 1; j < order.Len(); j++) {
		i = j;
		tempDat = order[j].Dat;
		while (i > 0 && order[i].Key == order[i - 1].Key && tempDat < order[i - 1].Dat) {
			order[i].Dat = order[i - 1].Dat;
			i--;
		}
		order[i].Dat = tempDat;
	}
}

void MultTempMotifCounter::FixSameTimestampOrder(TVec<TIntQuadruple>& order) {
	int i;
	TIntIntPrPr tempDat;
	for (int j = 1; j < order.Len(); j++) {
		i = j;
		tempDat = order[j].Dat;
		while (i > 0 && order[i].Key == order[i - 1].Key && tempDat.Val1 < order[i - 1].Dat.Val1) {
			order[i].Dat = order[i - 1].Dat;
			i--;
		}
		order[i].Dat = tempDat;
	}
}

void MultTempMotifCounter::Count3MTEdge2Node(int u, int v, double delta,
	Counter8D& counts) {

	// Sort event list by time
	TVec<TIntQuadruple> combined;
	int index = 0;
	AddStarEdges(combined, index, u, v, 0);
	AddStarEdges(combined, index, v, u, 1);
	combined.Sort();
	FixSameTimestampOrder(combined);

	// Get the counts
	TIntV in_out(combined.Len());
	TIntV timestamps(combined.Len());
	TIntV layers(combined.Len());

	//added attributes
	TInt attrSrc = attribute_data[u];
	TInt attrDst = attribute_data[v];

	for (int k = 0; k < combined.Len(); k++) {
		in_out[k] = combined[k].Dat.Val2.Val1;
		timestamps[k] = combined[k].Key;
		layers[k] = combined[k].Dat.Val2.Val2;
	}

	if (timestamps[0] == -1) {//partial_ ) {

		ThreeMpTEdgeMotifCounter pcounter(2, nrlayers_, nratts_);
		pcounter.Count(in_out, timestamps, layers, attrSrc, attrDst, delta, counts);
	}
	else {

		ThreeMTEdgeMotifCounter counter(2, nrlayers_, nratts_);
		counter.Count(in_out, timestamps, layers, attrSrc, attrDst, delta, counts);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Star counting methods
void MultTempMotifCounter::AddStarEdges(TVec<TIntQuadruple>& combined, int& index, int u, int v,
	int key) {
	if (HasEdges(u, v)) {
		const TIntV& timestamps = temporal_data_[u].GetDat(v);
		const TIntV& layers = layer_data_[u].GetDat(v);
		for (int i = 0; i < timestamps.Len(); i++) {
			combined.Add(TIntQuadruple(timestamps[i], TIntIntPrPr(index, TIntPr(key, layers[i]))));
			index++;
		}
	}
}

void MultTempMotifCounter::AddMultStarEdgeData(TVec<TIntPair>& ts_indices,
	TVec<MultStarEdgeData>& events,
	int& index, int u, int v, int nbr, int key, int att) {
	if (HasEdges(u, v)) {
		const TIntV& ts_vec = temporal_data_[u].GetDat(v);
		const TIntV& layers = layer_data_[u].GetDat(v);
		for (int j = 0; j < ts_vec.Len(); ++j) {
			ts_indices.Add(TIntPair(ts_vec[j], index));
			events.Add(MultStarEdgeData(nbr, key, layers[j], att));
			index++;
		}
	}
}

void MultTempMotifCounter::Count3MTEdge3NodeStars(double delta, Counter9D& pre_counts,
	Counter9D& pos_counts,
	Counter9D& mid_counts) {
	TIntV centers;
	GetAllNodes(centers);
	pre_counts = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	pos_counts = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	mid_counts = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	// Get counts for each node as the center
#pragma omp parallel for schedule(dynamic)  
	for (int c = 0; c < centers.Len(); c++) {
		// Gather all adjacent events
		int center = centers[c];
		TVec<TIntPair> ts_indices;
		TVec<MultStarEdgeData> events;
		//TNGraph::TNodeI NI = static_graph_->GetNI(center);
		int index = 0;
		TIntV nbrs;
		GetAllNeighbors(center, nbrs);
		int nbr_index = 0;
		for (; nbr_index < nbrs.Len(); nbr_index++) {
			int nbr = nbrs[nbr_index];
			AddMultStarEdgeData(ts_indices, events, index, center, nbr, nbr_index, 0, attribute_data[nbr]);
			AddMultStarEdgeData(ts_indices, events, index, nbr, center, nbr_index, 1, attribute_data[nbr]);
		}
		ts_indices.Sort();
		FixSameTimestampOrder(ts_indices);
		TIntV timestamps;
		TVec<MultStarEdgeData> ordered_events;
		for (int j = 0; j < ts_indices.Len(); j++) {
			timestamps.Add(ts_indices[j].Key);
			ordered_events.Add(events[ts_indices[j].Dat]);
		}

		if (timestamps[0] == -1) {
			ThreeMpTEdgeStarCounter ptesc(nbr_index + 1, nrlayers_, nratts_, attribute_data[center]);
			// dirs: outgoing --> 0, incoming --> 1
			ptesc.Count(ordered_events, timestamps, delta);
#pragma omp critical
			{ // Update counts
				for (int dir1 = 0; dir1 < 2; ++dir1) {
					for (int lay1 = 0; lay1 < nrlayers_; ++lay1) {
						for (int dir2 = 0; dir2 < 2; ++dir2) {
							for (int lay2 = 0; lay2 < nrlayers_; ++lay2) {
								for (int dir3 = 0; dir3 < 2; ++dir3) {
									for (int lay3 = 0; lay3 < nrlayers_; ++lay3) {
										for (int att1 = 0; att1 < nratts_; ++att1) {
											for (int att2 = 0; att2 < nratts_; ++att2) {
												for (int att3 = 0; att3 < nratts_; ++att3) {
													pre_counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3) +=
														ptesc.PreCount(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
													pos_counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3) +=
														ptesc.PosCount(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
													mid_counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3) +=
														ptesc.MidCount(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		else {
			ThreeMTEdgeStarCounter tesc(nbr_index + 1, nrlayers_, nratts_, attribute_data[center]);
			// dirs: outgoing --> 0, incoming --> 1
			tesc.Count(ordered_events, timestamps, delta);
#pragma omp critical
			{ // Update counts
				for (int dir1 = 0; dir1 < 2; ++dir1) {
					for (int lay1 = 0; lay1 < nrlayers_; ++lay1) {
						for (int dir2 = 0; dir2 < 2; ++dir2) {
							for (int lay2 = 0; lay2 < nrlayers_; ++lay2) {
								for (int dir3 = 0; dir3 < 2; ++dir3) {
									for (int lay3 = 0; lay3 < nrlayers_; ++lay3) {
										for (int att1 = 0; att1 < nratts_; ++att1) {
											for (int att2 = 0; att2 < nratts_; ++att2) {
												for (int att3 = 0; att3 < nratts_; ++att3) {
													pre_counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3) +=
														tesc.PreCount(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
													pos_counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3) +=
														tesc.PosCount(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
													mid_counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3) +=
														tesc.MidCount(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

		// Subtract off edge-wise counts
		for (int nbr_id = 0; nbr_id < nbrs.Len(); nbr_id++) {
			int nbr = nbrs[nbr_id];
			Counter8D edge_counts;
			Count3MTEdge2Node(center, nbr, delta, edge_counts);
#pragma omp critical
			{
				for (int dir1 = 0; dir1 < 2; ++dir1) {
					for (int lay1 = 0; lay1 < nrlayers_; ++lay1) {
						for (int dir2 = 0; dir2 < 2; ++dir2) {
							for (int lay2 = 0; lay2 < nrlayers_; ++lay2) {
								for (int lay3 = 0; lay3 < nrlayers_; ++lay3) {
									for (int att1 = 0; att1 < nratts_; ++att1) {
										for (int att2 = 0; att2 < nratts_; ++att2) {
											for (int att3 = 0; att3 < nratts_; ++att3) {
												pre_counts(dir1, lay1, dir2, lay2, 0, lay3, att1, att2, att2) -= edge_counts(dir1, lay1, dir2, lay2, 0, lay3, att1, att2);
												pre_counts(dir1, lay1, dir2, lay2, 1, lay3, att1, att2, att2) -= edge_counts(dir1, lay1, dir2, lay2, 1, lay3, att2, att1);
												pos_counts(dir1, lay1, dir2, lay2, 0, lay3, att1, att2, att2) -= edge_counts(dir1, lay1, dir2, lay2, 0, lay3, att1, att2);
												pos_counts(dir1, lay1, dir2, lay2, 1, lay3, att1, att2, att2) -= edge_counts(dir1, lay1, dir2, lay2, 1, lay3, att2, att1);
												mid_counts(dir1, lay1, dir2, lay2, 0, lay3, att1, att2, att2) -= edge_counts(dir1, lay1, dir2, lay2, 0, lay3, att1, att2);
												mid_counts(dir1, lay1, dir2, lay2, 1, lay3, att1, att2, att2) -= edge_counts(dir1, lay1, dir2, lay2, 1, lay3, att2, att1);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Triad counting methods
void MultTempMotifCounter::AddMultTriadEdgeData(TVec<MultTriadEdgeData>& events,
	TVec<TIntPair>& ts_indices,
	int& index, int u, int v, int nbr,
	int key1, int key2, int att) {
	if (HasEdges(u, v)) {
		const TIntV& timestamps = temporal_data_[u].GetDat(v);
		const TIntV& layers = layer_data_[u].GetDat(v);
		for (int i = 0; i < timestamps.Len(); i++) {
			ts_indices.Add(TIntPair(timestamps[i], index));
			events.Add(MultTriadEdgeData(nbr, key1, key2, layers[i], att));
			++index;
		}
	}
}

void MultTempMotifCounter::Count3MTEdgeTriads(double delta, Counter9D& counts) {
	counts = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);

	// Get the counts on each undirected edge
	TVec< THash<TInt, TInt> > edge_counts(static_graph_->GetMxNId());
	TVec< THash<TInt, TIntV> > assignments(static_graph_->GetMxNId());
	for (TNGraph::TEdgeI it = static_graph_->BegEI();
		it < static_graph_->EndEI(); it++) {
		int src = it.GetSrcNId();
		int dst = it.GetDstNId();
		int min_node = MIN(src, dst);
		int max_node = MAX(src, dst);
		edge_counts[min_node](max_node) += temporal_data_[src](dst).Len();
		assignments[min_node](max_node) = TIntV();
	}

	// Assign triangles to the edge with the most events
	TIntV Us, Vs, Ws;
	GetAllStaticTriangles(Us, Vs, Ws);
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Us.Len(); i++) {
		int u = Us[i];
		int v = Vs[i];
		int w = Ws[i];
		int counts_uv = edge_counts[MIN(u, v)].GetDat(MAX(u, v));
		int counts_uw = edge_counts[MIN(u, w)].GetDat(MAX(u, w));
		int counts_vw = edge_counts[MIN(v, w)].GetDat(MAX(v, w));
		if (counts_uv >= MAX(counts_uw, counts_vw)) {
#pragma omp critical
			{
				TIntV& assignment = assignments[MIN(u, v)].GetDat(MAX(u, v));
				assignment.Add(w);
			}
		}
		else if (counts_uw >= MAX(counts_uv, counts_vw)) {
#pragma omp critical
			{
				TIntV& assignment = assignments[MIN(u, w)].GetDat(MAX(u, w));
				assignment.Add(v);
			}
		}
		else if (counts_vw >= MAX(counts_uv, counts_uw)) {
#pragma omp critical
			{
				TIntV& assignment = assignments[MIN(v, w)].GetDat(MAX(v, w));
				assignment.Add(u);
			}
		}
	}

	TVec<TIntPair> all_edges;
	TIntV all_nodes;
	GetAllNodes(all_nodes);
	for (int node_id = 0; node_id < all_nodes.Len(); node_id++) {
		int u = all_nodes[node_id];
		TIntV nbrs;
		GetAllNeighbors(u, nbrs);
		for (int nbr_id = 0; nbr_id < nbrs.Len(); nbr_id++) {
			int v = nbrs[nbr_id];
			if (assignments[u].IsKey(v) && assignments[u].GetDat(v).Len() > 0) {
				all_edges.Add(TIntPair(u, v));
			}
		}
	}

	// Count triangles on edges with the assigned neighbours
#pragma omp parallel for schedule(dynamic)
	for (int edge_id = 0; edge_id < all_edges.Len(); edge_id++) {
		TIntPair edge = all_edges[edge_id];
		int u = edge.Key;
		int v = edge.Dat;
		// Continue if no assignment
		if (!assignments[u].IsKey(v)) { continue; }
		TIntV& uv_assignment = assignments[u].GetDat(v);
		// Continue if no data
		if (uv_assignment.Len() == 0) { continue; }
		// Get all events on (u, v)
		TVec<MultTriadEdgeData> events;
		TVec<TIntPair> ts_indices;
		int index = 0;
		int nbr_index = 0;
		// Assign indices from 0, 1, ..., num_nbrs + 2
		std::vector<int> nbr_attributes;
		AddMultTriadEdgeData(events, ts_indices, index, u, v, nbr_index, 0, 1, attribute_data[u]);
		nbr_index++;
		nbr_attributes.push_back(attribute_data[u]);
		AddMultTriadEdgeData(events, ts_indices, index, v, u, nbr_index, 0, 0, attribute_data[v]);
		nbr_attributes.push_back(attribute_data[v]);
		nbr_index++;
		// Get all events on triangles assigned to (u, v)
		for (int w_id = 0; w_id < uv_assignment.Len(); w_id++) {
			int w = uv_assignment[w_id];
			AddMultTriadEdgeData(events, ts_indices, index, w, u, nbr_index, 0, 0, attribute_data[w]);
			AddMultTriadEdgeData(events, ts_indices, index, w, v, nbr_index, 0, 1, attribute_data[w]);
			AddMultTriadEdgeData(events, ts_indices, index, u, w, nbr_index, 1, 0, attribute_data[w]);
			AddMultTriadEdgeData(events, ts_indices, index, v, w, nbr_index, 1, 1, attribute_data[w]);
			nbr_index++;
			nbr_attributes.push_back(attribute_data[w]);
		}
		// Put events in sorted order
		ts_indices.Sort();
		FixSameTimestampOrder(ts_indices);
		TIntV timestamps(ts_indices.Len());
		TVec<MultTriadEdgeData> sorted_events(ts_indices.Len());
		for (int i = 0; i < ts_indices.Len(); i++) {
			timestamps[i] = ts_indices[i].Key;
			sorted_events[i] = events[ts_indices[i].Dat];
		}

		// Get the counts and update the counter
		if (timestamps[0] == -1) {
			ThreeMpTEdgeTriadCounter ptetc(nbr_index, 0, 1, nrlayers_, nratts_, nbr_attributes, attribute_data[u], attribute_data[v]);
			ptetc.Count(sorted_events, timestamps, delta);
#pragma omp critical
			{
				for (int dir1 = 0; dir1 < 2; ++dir1) {
					for (int lay1 = 0; lay1 < nrlayers_; ++lay1) {
						for (int dir2 = 0; dir2 < 2; ++dir2) {
							for (int lay2 = 0; lay2 < nrlayers_; ++lay2) {
								for (int dir3 = 0; dir3 < 2; ++dir3) {
									for (int lay3 = 0; lay3 < nrlayers_; ++lay3) {
										for (int att1 = 0; att1 < nratts_; ++att1) {
											for (int att2 = 0; att2 < nratts_; ++att2) {
												for (int att3 = 0; att3 < nratts_; ++att3) {
													counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3) +=
														ptetc.Counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		else {
			ThreeMTEdgeTriadCounter tetc(nbr_index, 0, 1, nrlayers_, nratts_, nbr_attributes, attribute_data[u], attribute_data[v]);
			tetc.Count(sorted_events, timestamps, delta);
#pragma omp critical
			{
				for (int dir1 = 0; dir1 < 2; ++dir1) {
					for (int lay1 = 0; lay1 < nrlayers_; ++lay1) {
						for (int dir2 = 0; dir2 < 2; ++dir2) {
							for (int lay2 = 0; lay2 < nrlayers_; ++lay2) {
								for (int dir3 = 0; dir3 < 2; ++dir3) {
									for (int lay3 = 0; lay3 < nrlayers_; ++lay3) {
										for (int att1 = 0; att1 < nratts_; ++att1) {
											for (int att2 = 0; att2 < nratts_; ++att2) {
												for (int att3 = 0; att3 < nratts_; ++att3) {
													counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3) +=
														tetc.Counts(dir1, lay1, dir2, lay2, dir3, lay3, att1, att2, att3);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Generic three temporal edge motif counter
void ThreeMTEdgeMotifCounter::Count(const TIntV& event_string, const TIntV& timestamps,
	const TIntV& layers, const TInt& attSrc, const TInt& attDst, double delta, Counter8D& counts) {
	// Initialize everything to empty

	counts1_ = Counter2D(size_, nrlayers_);

	counts2_ = Counter4D(size_, nrlayers_, size_, nrlayers_);
	counts3_ = Counter8D(size_, nrlayers_, size_, nrlayers_, size_, nrlayers_, nratts_, nratts_);

	if (event_string.Len() != timestamps.Len() || event_string.Len() != layers.Len()) {
		TExcept::Throw("Number of events must equal number of timestamps and layers");
	}
	int start = 0;
	for (int end = 0; end < event_string.Len(); end++) {
		while (double(timestamps[start]) + delta < double(timestamps[end])) {
			DecrementCounts(event_string[start], layers[start]);
			start++;
		}
		IncrementCounts(event_string[end], layers[end], attSrc, attDst);
	}
	counts = counts3_;
}

void ThreeMTEdgeMotifCounter::IncrementCounts(int event, int layer, const TInt& attSrc, const TInt& attDst) {
	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			for (int k = 0; k < size_; k++) {
				for (int l = 0; l < nrlayers_; l++) {

					counts3_(i, j, k, l, event, layer, attSrc, attDst) += counts2_(i, j, k, l);//TODO update with extra loops
				}
			}
		}
	}
	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			counts2_(i, j, event, layer) += counts1_(i, j);
		}
	}
	counts1_(event, layer) += 1;
}

void ThreeMTEdgeMotifCounter::DecrementCounts(int event, int layer) {
	counts1_(event, layer)--;
	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < nrlayers_; j++)
			counts2_(event, layer, i, j) -= counts1_(i, j);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Generic three temporal edge motif counter for partial timing
void ThreeMpTEdgeMotifCounter::Count(const TIntV& event_string, const TIntV& timestamps,
	const TIntV& layers, const TInt& attrSrc, const TInt& attrDst, double delta, Counter8D& counts) {
	// Initialize everything to empty

	counts1_ = Counter2D(size_, nrlayers_);
	counts2_ = Counter4D(size_, nrlayers_, size_, nrlayers_);
	counts3_ = Counter8D(size_, nrlayers_, size_, nrlayers_, size_, nrlayers_, nratts_, nratts_);
	pcounts1_ = Counter2D(size_, nrlayers_);
	pcounts2_ = Counter4D(size_, nrlayers_, size_, nrlayers_);
	if (event_string.Len() != timestamps.Len() || event_string.Len() != layers.Len()) {
		TExcept::Throw("Number of events must equal number of timestamps and layers");
	}
	int end = 0;
	while (end < event_string.Len() && timestamps[end] == -1) {
		PreProcess(event_string[end], layers[end], attrSrc, attrDst); //TODO fix this
		end++;
	}
	int start = end;
	for (; end < event_string.Len(); end++) {
		while (double(timestamps[start]) + delta < double(timestamps[end])) {
			DecrementCounts(event_string[start], layers[start], attrSrc, attrDst);
			start++;
		}
		IncrementCounts(event_string[end], layers[end], attrSrc, attrDst);
	}
	counts = counts3_;
}

void ThreeMpTEdgeMotifCounter::PreProcess(int event, int layer, int att1, int att2) {
	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			for (int k = 0; k < size_; k++) {
				for (int l = 0; l < nrlayers_; l++)
					counts3_(i, j, k, l, event, layer, att1, att2) += pcounts2_(i, j, k, l);
			}
		}
	}
	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < nrlayers_; j++)
			pcounts2_(i, j, event, layer) += pcounts1_(i, j);
	}
	pcounts1_(event, layer) += 1;
}

void ThreeMpTEdgeMotifCounter::IncrementCounts(int event, int layer, int att1, int att2) {
	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			for (int k = 0; k < size_; k++) {
				for (int l = 0; l < nrlayers_; l++)
					counts3_(i, j, k, l, event, layer, att1, att2) += counts2_(i, j, k, l) + pcounts2_(i, j, k, l);
			}
		}
	}
	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < nrlayers_; j++)
			counts2_(i, j, event, layer) += counts1_(i, j) + pcounts1_(i, j);
	}
	counts1_(event, layer) += 1;
}

void ThreeMpTEdgeMotifCounter::DecrementCounts(int event, int layer, int att1, int att2) {
	counts1_(event, layer)--;
	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			counts2_(event, layer, i, j) -= counts1_(i, j);
			counts2_(i, j, event, layer) -= pcounts1_(i, j);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Generic three temporal edge, three node star and triad counter.
template <typename EdgeData>
void StarTriad3MTEdgeCounter<EdgeData>::Count(const TVec<EdgeData>& events,
	const TIntV& timestamps, double delta) {
	InitializeCounters();
	if (events.Len() != timestamps.Len()) {
		TExcept::Throw("Number of events must match number of timestamps.");
	}
	int start = 0;
	int end = 0;
	int L = timestamps.Len();
	for (int j = 0; j < L; j++) {
		double tj = double(timestamps[j]);
		// Adjust counts in pre-window [tj - delta, tj)
		while (start < L && double(timestamps[start]) < tj - delta) {
			PopPre(events[start]);
			start++;
		}
		// Adjust counts in post-window (tj, tj + delta]
		while (end < L && double(timestamps[end]) <= tj + delta) {
			PushPos(events[end]);
			end++;
		}
		// Move current event off post-window
		PopPos(events[j]);
		ProcessCurrent(events[j]);
		PushPre(events[j]);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Generic three temporal edge, three node star and triad counter for partial timing.
template <typename EdgeData>
void StarTriad3MpTEdgeCounter<EdgeData>::Count(const TVec<EdgeData>& events,
	const TIntV& timestamps, double delta) {
	InitializeCounters();
	if (events.Len() != timestamps.Len()) {
		TExcept::Throw("Number of events must match number of timestamps.");
	}
	int start = 0;
	int end = 0;
	int L = timestamps.Len();

	int j = 0;
	while (end < L) {//&& timestamps[end] == -1) {
		PreProcPushPos(events[end]);
		end++;
	}
	while (j < L && timestamps[j] == -1) {
		PreProcPopPos(events[j]);
		PreProcProcessCurrent(events[j]);
		PreProcPushPre(events[j]);
		j++;
	}
	start = j;
	end = j;
	for (; j < L; j++) {
		double tj = double(timestamps[j]);
		// Adjust counts in pre-window [tj - delta, tj)
		while (start < L && double(timestamps[start]) < tj - delta) {
			PopPre(events[start]);
			start++;
		}
		// Adjust counts in post-window (tj, tj + delta]
		while (end < L && double(timestamps[end]) <= tj + delta) {
			PushPos(events[end]);
			end++;
		}
		// Move current event off post-window
		PopPos(events[j]);
		ProcessCurrent(events[j]);
		PushPre(events[j]);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Methods for the main sub-routine in the fast star counting algorithm.
void ThreeMTEdgeStarCounter::InitializeCounters() {
	pre_sum_ = Counter5D(2, nrlayers_, 2, nrlayers_, nratts_);
	pos_sum_ = Counter5D(2, nrlayers_, 2, nrlayers_, nratts_);
	mid_sum_ = Counter5D(2, nrlayers_, 2, nrlayers_, nratts_);
	pre_counts_ = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	pos_counts_ = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	mid_counts_ = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	pre_nodes_ = Counter4D(2, nrlayers_, max_nodes_, nratts_);
	pos_nodes_ = Counter4D(2, nrlayers_, max_nodes_, nratts_);
}

void ThreeMTEdgeStarCounter::PopPre(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	pre_nodes_(dir, lay, nbr, att) -= 1;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			pre_sum_(dir, lay, i, j, att) -= pre_nodes_(i, j, nbr, att);
	}
}

void ThreeMTEdgeStarCounter::PopPos(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	pos_nodes_(dir, lay, nbr, att) -= 1;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			pos_sum_(dir, lay, i, j, att) -= pos_nodes_(i, j, nbr, att);
	}
}

void ThreeMTEdgeStarCounter::PushPre(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			pre_sum_(i, j, dir, lay, att) += pre_nodes_(i, j, nbr, att);
	}
	pre_nodes_(dir, lay, nbr, att) += 1;
}

void ThreeMTEdgeStarCounter::PushPos(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			pos_sum_(i, j, dir, lay, att) += pos_nodes_(i, j, nbr, att);
	}
	pos_nodes_(dir, lay, nbr, att) += 1;
}

void ThreeMTEdgeStarCounter::ProcessCurrent(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	// Decrement middle sum
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			mid_sum_(i, j, dir, lay, att) -= pre_nodes_(i, j, nbr, att);
		}
	}
	// Update counts
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < nrlayers_; l++) {
					for (int a = 0; a < nratts_; a++) {
						pre_counts_(i, j, k, l, dir, lay, c_att, att, a) += pre_sum_(i, j, k, l, a);
						pos_counts_(dir, lay, i, j, k, l, c_att, att, a) += pos_sum_(i, j, k, l, a);
						mid_counts_(i, j, dir, lay, k, l, c_att, att, a) += mid_sum_(i, j, k, l, a);

					}
				}
			}
		}
	}
	// Increment middle sum
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			mid_sum_(dir, lay, i, j, att) += pos_nodes_(i, j, nbr, att);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Methods for the main sub-routine in the fast partial timing star counting algorithm.
void ThreeMpTEdgeStarCounter::InitializeCounters() {
	pre_sum_ = Counter5D(2, nrlayers_, 2, nrlayers_, nratts_);
	pos_sum_ = Counter5D(2, nrlayers_, 2, nrlayers_, nratts_);
	mid_sum_ = Counter5D(2, nrlayers_, 2, nrlayers_, nratts_);
	pre_counts_ = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	pos_counts_ = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	mid_counts_ = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
	pre_nodes_ = Counter3D(2, nrlayers_, max_nodes_);
	pos_nodes_ = Counter3D(2, nrlayers_, max_nodes_);

	p_pre_mid_ = Counter4D(2, nrlayers_, 2, nrlayers_);

	p_pre_sum_ = Counter4D(2, nrlayers_, 2, nrlayers_);
	p_pos_sum_ = Counter4D(2, nrlayers_, 2, nrlayers_);
	p_mid_sum_ = Counter4D(2, nrlayers_, 2, nrlayers_);
	p_pre_nodes_ = Counter3D(2, nrlayers_, max_nodes_);
	p_pos_nodes_ = Counter3D(2, nrlayers_, max_nodes_);
}

void ThreeMpTEdgeStarCounter::PopPre(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	pre_nodes_(dir, lay, nbr) -= 1;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			pre_sum_(dir, lay, i, j, att) -= pre_nodes_(i, j, nbr);
			pre_sum_(i, j, dir, lay, att) -= p_pre_nodes_(i, j, nbr);
		}
	}
}

void ThreeMpTEdgeStarCounter::PopPos(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	pos_nodes_(dir, lay, nbr) -= 1;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			pos_sum_(dir, lay, i, j, att) -= pos_nodes_(i, j, nbr);
	}
}

void ThreeMpTEdgeStarCounter::PushPre(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			pre_sum_(i, j, dir, lay, att) += pre_nodes_(i, j, nbr) + p_pre_nodes_(i, j, nbr);
	}
	pre_nodes_(dir, lay, nbr) += 1;
}

void ThreeMpTEdgeStarCounter::PushPos(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			pos_sum_(i, j, dir, lay, att) += pos_nodes_(i, j, nbr);
			p_pre_mid_(i, j, dir, lay) += p_pre_nodes_(i, j, nbr);
		}
	}
	pos_nodes_(dir, lay, nbr) += 1;
}

void ThreeMpTEdgeStarCounter::ProcessCurrent(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	int att = event.att;
	// Decrement middle sum
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			mid_sum_(i, j, dir, lay, att) -= pre_nodes_(i, j, nbr);
			p_pre_mid_(i, j, dir, lay) -= p_pre_nodes_(i, j, nbr);
		}
	}
	// Update counts
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < nrlayers_; l++) {
					for (int a = 0; a < nratts_; a++) {
						pre_counts_(i, j, k, l, dir, lay, c_att, att, a) += pre_sum_(i, j, k, l, a);
						pos_counts_(dir, lay, i, j, k, l, c_att, att, a) += pos_sum_(i, j, k, l, a);
						mid_counts_(i, j, dir, lay, k, l, c_att, att, a) += mid_sum_(i, j, k, l, a);

					}
				}
			}
		}
	}
	// Increment middle sum
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			mid_sum_(dir, lay, i, j, att) += pos_nodes_(i, j, nbr);
	}
}

void ThreeMpTEdgeStarCounter::PreProcPopPos(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	p_pos_nodes_(dir, lay, nbr) -= 1;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			p_pos_sum_(dir, lay, i, j) -= p_pos_nodes_(i, j, nbr);
	}
}

void ThreeMpTEdgeStarCounter::PreProcPushPre(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			p_pre_sum_(i, j, dir, lay) += p_pre_nodes_(i, j, nbr);
	}
	p_pre_nodes_(dir, lay, nbr) += 1;
}

void ThreeMpTEdgeStarCounter::PreProcPushPos(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			p_pos_sum_(i, j, dir, lay) += p_pos_nodes_(i, j, nbr);
	}
	p_pos_nodes_(dir, lay, nbr) += 1;
}

void ThreeMpTEdgeStarCounter::PreProcProcessCurrent(const MultStarEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int lay = event.l;
	// Decrement middle sum
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			p_mid_sum_(i, j, dir, lay) -= p_pre_nodes_(i, j, nbr);
	}
	// Update counts
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++) {
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < nrlayers_; l++) {
				/*	pre_counts_(i, j, k, l, dir, lay) += p_pre_sum_(i, j, k, l);
					pos_counts_(dir, lay, i, j, k, l) += p_pos_sum_(i, j, k, l);
					mid_counts_(i, j, dir, lay, k, l) += p_mid_sum_(i, j, k, l);*/
				}
			}
		}
	}
	// Increment middle sum
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nrlayers_; j++)
			p_mid_sum_(dir, lay, i, j) += p_pos_nodes_(i, j, nbr);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Methods for the main sub-routine in the fast triangle counting algorithm.
void ThreeMTEdgeTriadCounter::InitializeCounters() {
	pre_nodes_ = Counter4D(2, nrlayers_, 2, max_nodes_);
	pos_nodes_ = Counter4D(2, nrlayers_, 2, max_nodes_);
	pre_sum_ = Counter6D(2, 2, nrlayers_, 2, nrlayers_, nratts_);
	pos_sum_ = Counter6D(2, 2, nrlayers_, 2, nrlayers_, nratts_);
	mid_sum_ = Counter6D(2, 2, nrlayers_, 2, nrlayers_, nratts_);
	triad_counts_ = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);
}

void ThreeMTEdgeTriadCounter::PopPre(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	if (!IsEdgeNode(nbr)) {
		pre_nodes_(dir, lay, u_or_v, nbr) -= 1;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				pre_sum_(u_or_v, dir, lay, i, j, att) -= pre_nodes_(i, j, 1 - u_or_v, nbr);
		}
	}
}

void ThreeMTEdgeTriadCounter::PopPos(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	if (!IsEdgeNode(nbr)) {
		pos_nodes_(dir, lay, u_or_v, nbr) -= 1;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				pos_sum_(u_or_v, dir, lay, i, j, att) -= pos_nodes_(i, j, 1 - u_or_v, nbr);
		}
	}
}

void ThreeMTEdgeTriadCounter::PushPre(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				pre_sum_(1 - u_or_v, i, j, dir, lay, att) += pre_nodes_(i, j, 1 - u_or_v, nbr);
		}
		pre_nodes_(dir, lay, u_or_v, nbr) += 1;
	}
}

void ThreeMTEdgeTriadCounter::PushPos(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				pos_sum_(1 - u_or_v, i, j, dir, lay, att) += pos_nodes_(i, j, 1 - u_or_v, nbr);
		}
		pos_nodes_(dir, lay, u_or_v, nbr) += 1;
	}
}

void ThreeMTEdgeTriadCounter::ProcessCurrent(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l; 
	int att = attribute_data[nbr];
	// Adjust middle sums
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++) {
				mid_sum_(1 - u_or_v, i, j, dir, lay, att) -= pre_nodes_(i, j, 1 - u_or_v, nbr);
				mid_sum_(u_or_v, dir, lay, i, j, att) += pos_nodes_(i, j, 1 - u_or_v, nbr);
			}
		}
	}
	// Update counts
	if (IsEdgeNode(nbr)) {
		// Determine if the event edge is u --> v or v --> u
		int u_to_v = 0;
		if (((nbr == node_u_) && dir == 0) || ((nbr == node_v_) && dir == 1)) {
			u_to_v = 1;
			int temp = u_att;
			u_att = v_att;
			v_att = temp;
		}

		for (int i = 0; i < nrlayers_; i++) {
			for (int j = 0; j < nrlayers_; j++) {
				for (int k = 0; k < nratts_; k++) {
					// i --> j, k --> j, i --> k    
					triad_counts_(0, i, 0, lay, 0, j, k, u_att, v_att) += mid_sum_(u_to_v, 0, i, 0, j, k);
					triad_counts_(0, lay, 0, i, 0, j, u_att, v_att, k) += pos_sum_(u_to_v, 0, i, 1, j, k);
					triad_counts_(0, i, 0, j, 0, lay, u_att, k, v_att) += pre_sum_(1 - u_to_v, 1, i, 1, j, k);
					// i --> j, k --> i, j --> k
					triad_counts_(1, i, 0, lay, 0, j, k, u_att, v_att) += mid_sum_(u_to_v, 1, i, 0, j, k);
					triad_counts_(1, lay, 0, i, 0, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 0, i, 1, j, k);
					triad_counts_(1, i, 0, j, 0, lay, v_att, k, u_att) += pre_sum_(1 - u_to_v, 0, i, 1, j, k);
					// i --> j, j --> k, i --> k
					triad_counts_(0, i, 1, lay, 0, j, k, v_att, u_att) += mid_sum_(1 - u_to_v, 0, i, 0, j, k);
					triad_counts_(0, lay, 1, i, 0, j, u_att, v_att, k) += pos_sum_(u_to_v, 1, i, 1, j, k);
					triad_counts_(0, i, 1, j, 0, lay, u_att, k, v_att) += pre_sum_(1 - u_to_v, 1, i, 0, j, k);
					// i --> j, i --> k, j --> k
					triad_counts_(1, i, 1, lay, 0, j, k, v_att, u_att) += mid_sum_(1 - u_to_v, 1, i, 0, j, k);
					triad_counts_(1, lay, 1, i, 0, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 1, i, 1, j, k);
					triad_counts_(1, i, 1, j, 0, lay, v_att, k, u_att) += pre_sum_(1 - u_to_v, 0, i, 0, j, k);
					// i --> j, k --> j, k --> i
					triad_counts_(0, i, 0, lay, 1, j, u_att, k, v_att) += mid_sum_(u_to_v, 0, i, 1, j, k);
					triad_counts_(0, lay, 0, i, 1, j, u_att, v_att, k) += pos_sum_(u_to_v, 0, i, 0, j, k);
					triad_counts_(0, i, 0, j, 1, lay, k, u_att, v_att) += pre_sum_(u_to_v, 1, i, 1, j, k);
					// i --> j, k --> i, k --> j
					triad_counts_(1, i, 0, lay, 1, j, u_att, k, v_att) += mid_sum_(u_to_v, 1, i, 1, j, k);
					triad_counts_(1, lay, 0, i, 1, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 0, i, 0, j, k);
					triad_counts_(1, i, 0, j, 1, lay, k, v_att, u_att) += pre_sum_(u_to_v, 0, i, 1, j, k);
					// i --> j, j --> k, k --> i
					triad_counts_(0, i, 1, lay, 1, j, v_att, k, u_att) += mid_sum_(1 - u_to_v, 0, i, 1, j, k);
					triad_counts_(0, lay, 1, i, 1, j, u_att, v_att, k) += pos_sum_(u_to_v, 1, i, 0, j, k);
					triad_counts_(0, i, 1, j, 1, lay, k, u_att, v_att) += pre_sum_(u_to_v, 1, i, 0, j, k);
					// i --> j, i --> k, k --> j
					triad_counts_(1, i, 1, lay, 1, j, v_att, k, u_att) += mid_sum_(1 - u_to_v, 1, i, 1, j, k);
					triad_counts_(1, lay, 1, i, 1, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 1, i, 0, j, k);
					triad_counts_(1, i, 1, j, 1, lay, k, v_att, u_att) += pre_sum_(u_to_v, 0, i, 0, j, k);
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Methods for the main sub-routine in the fast partial timing triangle counting algorithm.
void ThreeMpTEdgeTriadCounter::InitializeCounters() {
	pre_nodes_ = Counter4D(2, nrlayers_, 2, max_nodes_);
	pos_nodes_ = Counter4D(2, nrlayers_, 2, max_nodes_);
	pre_sum_ = Counter6D(2, 2, nrlayers_, 2, nrlayers_, nratts_);
	pos_sum_ = Counter6D(2, 2, nrlayers_, 2, nrlayers_, nratts_);
	mid_sum_ = Counter6D(2, 2, nrlayers_, 2, nrlayers_, nratts_);
	triad_counts_ = Counter9D(2, nrlayers_, 2, nrlayers_, 2, nrlayers_, nratts_, nratts_, nratts_);

	p_pre_mid_ = Counter5D(2, 2, nrlayers_, 2, nrlayers_);

	p_pre_nodes_ = Counter4D(2, nrlayers_, 2, max_nodes_);
	p_pos_nodes_ = Counter4D(2, nrlayers_, 2, max_nodes_);
	p_pre_sum_ = Counter5D(2, 2, nrlayers_, 2, nrlayers_);
	p_pos_sum_ = Counter5D(2, 2, nrlayers_, 2, nrlayers_);
	p_mid_sum_ = Counter5D(2, 2, nrlayers_, 2, nrlayers_);
}

void ThreeMpTEdgeTriadCounter::PopPre(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	if (!IsEdgeNode(nbr)) {
		pre_nodes_(dir, lay, u_or_v, nbr) -= 1;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++) {
				pre_sum_(u_or_v, dir, lay, i, j, att) -= pre_nodes_(i, j, 1 - u_or_v, nbr);
				pre_sum_(1 - u_or_v, i, j, dir, lay, att) -= p_pre_nodes_(i, j, 1 - u_or_v, nbr);
			}
		}
	}
}

void ThreeMpTEdgeTriadCounter::PopPos(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	if (!IsEdgeNode(nbr)) {
		pos_nodes_(dir, lay, u_or_v, nbr) -= 1;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				pos_sum_(u_or_v, dir, lay, i, j, att) -= pos_nodes_(i, j, 1 - u_or_v, nbr);
		}
	}
}

void ThreeMpTEdgeTriadCounter::PushPre(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				pre_sum_(1 - u_or_v, i, j, dir, lay, att) += pre_nodes_(i, j, 1 - u_or_v, nbr)
				+ p_pre_nodes_(i, j, 1 - u_or_v, nbr);
		}
		pre_nodes_(dir, lay, u_or_v, nbr) += 1;
	}
}

void ThreeMpTEdgeTriadCounter::PushPos(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++) {
				pos_sum_(1 - u_or_v, i, j, dir, lay, att) += pos_nodes_(i, j, 1 - u_or_v, nbr);
				p_pre_mid_(1 - u_or_v, i, j, dir, lay) += p_pre_nodes_(i, j, 1 - u_or_v, nbr);
			}
		}
		pos_nodes_(dir, lay, u_or_v, nbr) += 1;
	}
}

void ThreeMpTEdgeTriadCounter::ProcessCurrent(const MultTriadEdgeData& event) {

	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	int att = attribute_data[nbr];
	// Adjust middle sums
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++) {
				mid_sum_(1 - u_or_v, i, j, dir, lay, att) -= pre_nodes_(i, j, 1 - u_or_v, nbr);
				mid_sum_(u_or_v, dir, lay, i, j, att) += pos_nodes_(i, j, 1 - u_or_v, nbr);
			}
		}
	}
	// Update counts
	if (IsEdgeNode(nbr)) {
		// Determine if the event edge is u --> v or v --> u
		int u_to_v = 0;
		if (((nbr == node_u_) && dir == 0) || ((nbr == node_v_) && dir == 1)) {
			u_to_v = 1;
		}
		if (u_to_v == 1) {
			int temp = u_att;
			u_att = v_att;
			v_att = temp;
		}

		for (int i = 0; i < nrlayers_; i++) {
			for (int j = 0; j < nrlayers_; j++) {
				for (int k = 0; k < nratts_; k++) {
					// i --> j, k --> j, i --> k    
					triad_counts_(0, i, 0, lay, 0, j, u_att, v_att, k) += mid_sum_(u_to_v, 0, i, 0, j, k)
						+ p_pre_mid_(u_to_v, 0, i, 0, j);
					triad_counts_(0, lay, 0, i, 0, j, u_att, v_att, k) += pos_sum_(u_to_v, 0, i, 1, j, k);
					triad_counts_(0, i, 0, j, 0, lay, u_att, v_att, k) += pre_sum_(1 - u_to_v, 1, i, 1, j, k)
						+ p_pre_sum_(1 - u_to_v, 1, i, 1, j);
					// i --> j, k --> i, j --> k
					triad_counts_(1, i, 0, lay, 0, j, u_att, v_att, k) += mid_sum_(u_to_v, 1, i, 0, j, k)
						+ p_pre_mid_(u_to_v, 1, i, 0, j);
					triad_counts_(1, lay, 0, i, 0, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 0, i, 1, j, k);
					triad_counts_(1, i, 0, j, 0, lay, u_att, v_att, k) += pre_sum_(1 - u_to_v, 0, i, 1, j, k)
						+ p_pre_sum_(1 - u_to_v, 0, i, 1, j);
					// i --> j, j --> k, i --> k
					triad_counts_(0, i, 1, lay, 0, j, u_att, v_att, k) += mid_sum_(1 - u_to_v, 0, i, 0, j, k)
						+ p_pre_mid_(1 - u_to_v, 0, i, 0, j);
					triad_counts_(0, lay, 1, i, 0, j, u_att, v_att, k) += pos_sum_(u_to_v, 1, i, 1, j, k);
					triad_counts_(0, i, 1, j, 0, lay, u_att, v_att, k) += pre_sum_(1 - u_to_v, 1, i, 0, j, k)
						+ p_pre_sum_(1 - u_to_v, 1, i, 0, j);
					// i --> j, i --> k, j --> k
					triad_counts_(1, i, 1, lay, 0, j, u_att, v_att, k) += mid_sum_(1 - u_to_v, 1, i, 0, j, k)
						+ p_pre_mid_(1 - u_to_v, 1, i, 0, j);
					triad_counts_(1, lay, 1, i, 0, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 1, i, 1, j, k);
					triad_counts_(1, i, 1, j, 0, lay, u_att, v_att, k) += pre_sum_(1 - u_to_v, 0, i, 0, j, k)
						+ p_pre_sum_(1 - u_to_v, 0, i, 0, j);
					// i --> j, k --> j, k --> i
					triad_counts_(0, i, 0, lay, 1, j, u_att, v_att, k) += mid_sum_(u_to_v, 0, i, 1, j, k)
						+ p_pre_mid_(u_to_v, 0, i, 1, j);
					triad_counts_(0, lay, 0, i, 1, j, u_att, v_att, k) += pos_sum_(u_to_v, 0, i, 0, j, k);
					triad_counts_(0, i, 0, j, 1, lay, u_att, v_att, k) += pre_sum_(u_to_v, 1, i, 1, j, k)
						+ p_pre_sum_(u_to_v, 1, i, 1, j);
					// i --> j, k --> i, k --> j
					triad_counts_(1, i, 0, lay, 1, j, u_att, v_att, k) += mid_sum_(u_to_v, 1, i, 1, j, k)
						+ p_pre_mid_(u_to_v, 1, i, 1, j);
					triad_counts_(1, lay, 0, i, 1, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 0, i, 0, j, k);
					triad_counts_(1, i, 0, j, 1, lay, u_att, v_att, k) += pre_sum_(u_to_v, 0, i, 1, j, k)
						+ p_pre_sum_(u_to_v, 0, i, 1, j);
					// i --> j, j --> k, k --> i
					triad_counts_(0, i, 1, lay, 1, j, u_att, v_att, k) += mid_sum_(1 - u_to_v, 0, i, 1, j, k)
						+ p_pre_mid_(1 - u_to_v, 0, i, 1, j);
					triad_counts_(0, lay, 1, i, 1, j, u_att, v_att, k) += pos_sum_(u_to_v, 1, i, 0, j, k);
					triad_counts_(0, i, 1, j, 1, lay, u_att, v_att, k) += pre_sum_(u_to_v, 1, i, 0, j, k)
						+ p_pre_sum_(u_to_v, 1, i, 0, j);
					// i --> j, i --> k, k --> j
					triad_counts_(1, i, 1, lay, 1, j, u_att, v_att, k) += mid_sum_(1 - u_to_v, 1, i, 1, j, k)
						+ p_pre_mid_(1 - u_to_v, 1, i, 1, j);
					triad_counts_(1, lay, 1, i, 1, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 1, i, 0, j, k);
					triad_counts_(1, i, 1, j, 1, lay, u_att, v_att, k) += pre_sum_(u_to_v, 0, i, 0, j, k)
						+ p_pre_sum_(u_to_v, 0, i, 0, j);
				}
			}
		}
	}
}

void ThreeMpTEdgeTriadCounter::PreProcPopPos(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	if (!IsEdgeNode(nbr)) {
		p_pos_nodes_(dir, lay, u_or_v, nbr) -= 1;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				p_pos_sum_(u_or_v, dir, lay, i, j) -= p_pos_nodes_(i, j, 1 - u_or_v, nbr);
		}
	}
}

void ThreeMpTEdgeTriadCounter::PreProcPushPre(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				p_pre_sum_(1 - u_or_v, i, j, dir, lay) += p_pre_nodes_(i, j, 1 - u_or_v, nbr);
		}
		p_pre_nodes_(dir, lay, u_or_v, nbr) += 1;
	}
}

void ThreeMpTEdgeTriadCounter::PreProcPushPos(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++)
				p_pos_sum_(1 - u_or_v, i, j, dir, lay) += p_pos_nodes_(i, j, 1 - u_or_v, nbr);
		}
		p_pos_nodes_(dir, lay, u_or_v, nbr) += 1;
	}
}

void ThreeMpTEdgeTriadCounter::PreProcProcessCurrent(const MultTriadEdgeData& event) {
	int nbr = event.nbr;
	int dir = event.dir;
	int u_or_v = event.u_or_v;
	int lay = event.l;
	// Adjust middle sums
	if (!IsEdgeNode(nbr)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nrlayers_; j++) {
				p_mid_sum_(1 - u_or_v, i, j, dir, lay) -= p_pre_nodes_(i, j, 1 - u_or_v, nbr);
				p_mid_sum_(u_or_v, dir, lay, i, j) += p_pos_nodes_(i, j, 1 - u_or_v, nbr);
			}
		}
	}
	// Update counts
	if (IsEdgeNode(nbr)) {
		// Determine if the event edge is u --> v or v --> u
		int u_to_v = 0;
		if (((nbr == node_u_) && dir == 0) || ((nbr == node_v_) && dir == 1)) {
			u_to_v = 1;
		}

		int u_att = 0;
		int v_att = 0;
		int nbr_att = attribute_data[nbr];

		for (int i = 0; i < nrlayers_; i++) {
			for (int j = 0; j < nrlayers_; j++) {
				for (int k = 0; k < nratts_; k++) {
					// i --> j, k --> j, i --> k    
					triad_counts_(0, i, 0, lay, 0, j, k, u_att, v_att) += mid_sum_(u_to_v, 0, i, 0, j, k);
					triad_counts_(0, lay, 0, i, 0, j, u_att, v_att, k) += pos_sum_(u_to_v, 0, i, 1, j, k);
					triad_counts_(0, i, 0, j, 0, lay, u_att, k, v_att) += pre_sum_(1 - u_to_v, 1, i, 1, j, k);
					// i --> j, k --> i, j --> k
					triad_counts_(1, i, 0, lay, 0, j, k, u_att, v_att) += mid_sum_(u_to_v, 1, i, 0, j, k);
					triad_counts_(1, lay, 0, i, 0, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 0, i, 1, j, k);
					triad_counts_(1, i, 0, j, 0, lay, v_att, k, u_att) += pre_sum_(1 - u_to_v, 0, i, 1, j, k);
					// i --> j, j --> k, i --> k
					triad_counts_(0, i, 1, lay, 0, j, k, v_att, u_att) += mid_sum_(1 - u_to_v, 0, i, 0, j, k);
					triad_counts_(0, lay, 1, i, 0, j, u_att, v_att, k) += pos_sum_(u_to_v, 1, i, 1, j, k);
					triad_counts_(0, i, 1, j, 0, lay, u_att, k, v_att) += pre_sum_(1 - u_to_v, 1, i, 0, j, k);
					// i --> j, i --> k, j --> k
					triad_counts_(1, i, 1, lay, 0, j, k, v_att, u_att) += mid_sum_(1 - u_to_v, 1, i, 0, j, k);
					triad_counts_(1, lay, 1, i, 0, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 1, i, 1, j, k);
					triad_counts_(1, i, 1, j, 0, lay, v_att, k, u_att) += pre_sum_(1 - u_to_v, 0, i, 0, j, k);
					// i --> j, k --> j, k --> i
					triad_counts_(0, i, 0, lay, 1, j, u_att, k, v_att) += mid_sum_(u_to_v, 0, i, 1, j, k);
					triad_counts_(0, lay, 0, i, 1, j, u_att, v_att, k) += pos_sum_(u_to_v, 0, i, 0, j, k);
					triad_counts_(0, i, 0, j, 1, lay, k, u_att, v_att) += pre_sum_(u_to_v, 1, i, 1, j, k);
					// i --> j, k --> i, k --> j
					triad_counts_(1, i, 0, lay, 1, j, u_att, k, v_att) += mid_sum_(u_to_v, 1, i, 1, j, k);
					triad_counts_(1, lay, 0, i, 1, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 0, i, 0, j, k);
					triad_counts_(1, i, 0, j, 1, lay, k, v_att, u_att) += pre_sum_(u_to_v, 0, i, 1, j, k);
					// i --> j, j --> k, k --> i
					triad_counts_(0, i, 1, lay, 1, j, v_att, k, u_att) += mid_sum_(1 - u_to_v, 0, i, 1, j, k);
					triad_counts_(0, lay, 1, i, 1, j, u_att, v_att, k) += pos_sum_(u_to_v, 1, i, 0, j, k);
					triad_counts_(0, i, 1, j, 1, lay, k, u_att, v_att) += pre_sum_(u_to_v, 1, i, 0, j, k);
					// i --> j, i --> k, k --> j
					triad_counts_(1, i, 1, lay, 1, j, v_att, k, u_att) += mid_sum_(1 - u_to_v, 1, i, 1, j, k);
					triad_counts_(1, lay, 1, i, 1, j, u_att, v_att, k) += pos_sum_(1 - u_to_v, 1, i, 0, j, k);
					triad_counts_(1, i, 1, j, 1, lay, k, v_att, u_att) += pre_sum_(u_to_v, 0, i, 0, j, k);
				}
			}
		}
	}
}
