#include <cmath>
#include <random>
#include "low_disc_seq.h"

// this code is from: https://observablehq.com/@techsparx/an-improvement-on-bridsons-algorithm-for-poisson-disc-samp/2
bool poisson_far(poisson_grid_t& pg, float x, float y, float cell_size, float rad2) {
	float h = pg.size();
	float w = pg[0].size();
	int i = x / cell_size;
	int j = y / cell_size;
	int i0 = i - 2 > 0 ? i - 2 : 0;
	int j0 = j - 2 > 0 ? j - 2 : 0;
	int i1 = i + 3 < w ? i + 3 : w;
	int j1 = j + 3 < h ? j + 3 : h;

	for (j = j0; j < j1; ++j) {
		for (i = i0; i < i1; ++i) {
			if (pg[j][i].first != -1) {
				float dx = pg[j][i].first - x;
				float dy = pg[j][i].second - y;
				if (dx * dx + dy * dy < rad2) {
					return false;
				}
			}
		}
	}
	return true;
}
poisson_grid_t poission_disc_sampling(float width, float height, float radius) {
	int k = 4;
	const float rad_inner2 = radius * radius;
	float cell_size = radius / sqrt(2);
	int grid_width = ceil(width / cell_size);
	int grid_height = ceil(height / cell_size);

	poisson_grid_t pg(grid_height);
	for (auto& u : pg) {
		u = poisson_row_t(grid_width, { -1, -1 });
	}
	std::pair<float, float> s;
	s = { width / 2, height / 2 };
	poisson_row_t queue;
	pg[s.second / cell_size][s.first / cell_size] = s;
	queue.push_back(s);

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<float> urd_r(radius, 2 * radius);
	std::uniform_real_distribution<float> urd_a(0, 1);

	while (!queue.empty()) {
		int i = rand() % queue.size();
		s = queue[i];
		bool found = false;
		
		for (int j = 0; j < k; ++j) {
			float b1 = (rand() % 1000) / 1000.f;
			float b2 = (rand() % 1000) / 1000.f;
			float a = 2 * 3.1415926536f * b1;
			float r = radius + radius * b2;
			float x = s.first  + r * cos(a);
			float y = s.second + r * sin(a);
			if (0 <= x && x < width && 0 <= y && y < height && poisson_far(pg, x, y, cell_size, rad_inner2)) {
				found = true;
				pg[y / cell_size][x / cell_size] = {x, y};
				queue.push_back({ x, y });
				break;
			}
		}

		if (!found) {
			queue[i] = queue.back();
			queue.erase(queue.end() - 1);
		}
	}
	return pg;
}

