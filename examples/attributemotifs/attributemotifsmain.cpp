// attributemotifsmain.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "attributemotifs.h"
#include <time.h>
#include <iostream>
#include <string>

#ifdef USE_OPENMP
#include <omp.h>
#endif

struct timespec start, finish, startt, finishh;
double elapsed;

void attributeCounting(const TStr multilayer_temporal_graph_filename, const TStr output, const TFlt delta, const bool single_sim, const bool single_att) {

	MultTempMotifCounter mtmc(multilayer_temporal_graph_filename, single_sim, single_att);
	int total_motifs = (mtmc.nrlayers_ * mtmc.nrlayers_ * mtmc.nrlayers_);
	int total_attributes = (mtmc.nratts_ * mtmc.nratts_ * mtmc.nratts_);
	Counter3D counts2 = Counter3D(4, total_motifs, mtmc.nratts_ * mtmc.nratts_);
	Counter4D counts3 = Counter4D(4, 8, total_motifs, total_attributes);
	mtmc.Count3MTEdge23Node(delta, counts2, counts3);

	FILE* output_file = fopen(output.CStr(), "wt");
	fprintf(output_file, "2-Node-3-Edge Motifs:\n\n");

	for (int l = 0; l < mtmc.nratts_ * mtmc.nratts_; l++) {
		int q = (l / mtmc.nratts_);
		int r = l % mtmc.nratts_;
		fprintf(output_file, "Attribute permutation %d (", l + 1);
		fprintf(output_file, "%d, ", q);
		fprintf(output_file, "%d):\n", r);
		for (int k = 0; k < total_motifs; k++){
			int c = (k / (total_motifs * total_motifs));
			int b = (k / total_motifs) % total_motifs;
			int a = k % total_motifs;
			fprintf(output_file, "Layer permutation %d (", k + 1);
			fprintf(output_file, "%d, ", a);
			fprintf(output_file, "%d, ", b);
			fprintf(output_file, "%d):\n", c);
			for (int i = 0; i < 4; i++) {
				TUInt64 count = counts2(i, k, l);
				fprintf(output_file, "%s", (count.GetStr()).CStr());
				if (i < 3) { fprintf(output_file, " "); }
			}
			fprintf(output_file, "\n\n");
		}
	}

	fprintf(output_file, "3-Node-3-Edge Motifs:\n\n");

	for (int l = 0; l < total_attributes; l++) {
		int p = (l / (mtmc.nratts_ * mtmc.nratts_));
		int q = (l / mtmc.nratts_) % mtmc.nratts_;
		int r = l % mtmc.nratts_;
		fprintf(output_file, "Attribute permutation %d (", l + 1);
		fprintf(output_file, "%d, ", p);
		fprintf(output_file, "%d, ", q);
		fprintf(output_file, "%d):\n", r);
		for (int k = 0; k < total_motifs; k++) {
			int c = (k / (total_motifs * total_motifs));
			int b = (k / total_motifs) % total_motifs;
			int a = k % total_motifs;
			fprintf(output_file, "Layer permutation %d (", k + 1);
			fprintf(output_file, "%d, ", a);
			fprintf(output_file, "%d, ", b);
			fprintf(output_file, "%d):\n", c);
			for (int j = 0; j < 8; j++) {
				for (int i = 0; i < 4; i++) {
					TUInt64 count = counts3(i, j, k, l);
					fprintf(output_file, "%s", (count.GetStr()).CStr());
					if (i < 3) { fprintf(output_file, " "); }
				}
				fprintf(output_file, "\n");
				if (j == 1) {
					fprintf(output_file, "\n");

				}
			}
			fprintf(output_file, "\n\n");
		}
	}
}

int main(int argc, char* argv[]) {
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("attributemotifs. build: %s, %s. Time: %s",
		__TIME__, __DATE__, TExeTm::GetCurTm()));
	

		const TStr graph_filename =
		Env.GetIfArgPrefixStr("-i:", "example-annotated-temporal-graph.txt",
			"Input file base");
	const TStr output_filename =
		Env.GetIfArgPrefixStr("-o:", "example-output.txt",
				"Output file base in which to write counts");
	const TFlt delta =
		Env.GetIfArgPrefixFlt("-delta:", 4096, "Time window delta");
	const int num_threads =
		Env.GetIfArgPrefixInt("-nt:", 4, "Number of threads for parallelization");
	const bool single_sim =
		Env.GetIfArgPrefixBool("-sim:", false, "Whether to simulate all layers as a single layer (F or T)");
	const bool single_att =
		Env.GetIfArgPrefixBool("-att:", false, "Whether to simulate all attributes as a single attribute (F or T)");


#ifdef USE_OPENMP
	omp_set_num_threads(num_threads);
#endif


	attributeCounting(graph_filename, output_filename, delta, single_sim, single_att);

	return 0;
}
