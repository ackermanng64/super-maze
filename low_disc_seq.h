#pragma once
#include <vector>
// low-discrepency-sequence

using poisson_row_t = std::vector<std::pair<float, float>>;
using poisson_grid_t = std::vector <poisson_row_t>;
bool poisson_far(poisson_grid_t& v, float x, float y, float cell_size, float rad2);
poisson_grid_t poission_disc_sampling(float width, float height, float radius);