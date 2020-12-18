#include "geometry.h"
#include "voronoi.h"

std::vector<std::vector<int>> voronoi(const std::vector<tri>& tris) {
	std::vector<std::vector<int>> out;
	for (int i = 0; i < tris.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			bool shares_tri = false;
			for (int k = 0; k < tris.size(); ++k) {
				if (k == i) continue;
				for (int l = 0; l < 3; ++l) {
					bool b1 = (tris[i].a[j] == tris[k].a[l] && tris[i].a[(j + 1) % 3] == tris[k].a[(l + 1) % 3]);
					bool b2 = (tris[i].a[(j + 1) % 3] == tris[k].a[l] && tris[i].a[j] == tris[k].a[(l + 1) % 3]);
					if (b1 || b2) {
						out.push_back({ i, j, k});
						shares_tri = true;
						break;
					}
				}
				if (shares_tri) {
					break;
				}
			}
			if (!shares_tri) {
				out.push_back({ i, j, -1});
			}
		}
	}
	return out;
}