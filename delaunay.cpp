#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <stack>
#include "geometry.h"
#include "delaunay.h"


std::vector<tri> bowyer_watson_triangualtion(float* points, int N) {
	float mmx[2];
	float mmy[2];
	get_min_max_pos(mmx, mmy, points, N, 10);
	float a = 5 * (mmx[1] - mmx[0] > mmy[1] - mmy[0] ? mmy[1] - mmy[0] : mmx[1] - mmx[0]);
	float p[3][2] = { {mmx[0] - a, mmy[0] - a}, { mmx[0] + 2 * a, mmy[0] - a},  { mmx[0] + a / 2.0, mmy[0] + a * 1.80277564f} };
	float* s[3];
	
	std::vector<tri> triangles;
	std::vector<tri> gud_triangles;
	std::vector<int> bad_triangles;
	std::vector<std::pair<int, int>> polygon;

	tri t;
	t.a[0] = -3;
	t.a[1] = -2;
	t.a[2] = -1;
	
	triangles.push_back(t);
	
	for (int i = 0; i < N; ++i) {
		bad_triangles.clear();
		gud_triangles.clear();
		for (int j = 0; j < triangles.size(); ++j) {
			auto& u = triangles[j];
			float q[2];
			// mid
			for (int k = 0; k < 3; ++k) {
				if (u.a[k] < 0) {
					s[k] = p[u.a[k] + 3];
				}
				else {
					s[k] = &points[2 * u.a[k]];
				}
			}
			float c[2];
			get_tri_circumcirlce(s[0], s[1], s[2], c);
			q[0] = s[0][0] - c[0];
			q[1] = s[0][1] - c[1];
			c[0] = c[0] - points[2 * i];
			c[1] = c[1] - points[2 * i + 1];
			if (sqr_len(q) > sqr_len(c)) {
				bad_triangles.push_back(j);
				//printf("bad1: (%d, %d, %d)\n", triangles[j].a[0], triangles[j].a[1], triangles[j].a[2]);
			}
			else {
				gud_triangles.push_back(triangles[j]);
			}
		}
		polygon.clear();
		for (int j = 0; j < bad_triangles.size(); ++j) {
			auto& u = triangles[bad_triangles[j]];
			for (int k = 0; k < 3; ++k) {
				bool bad_edge = false;
				for (int l = 0; l < bad_triangles.size(); ++l) {
					if (l != j) {
						auto& v = triangles[bad_triangles[l]];
						for (int m = 0; m < 3; ++m) {
							if ((u.a[k] == v.a[m] && u.a[(k + 1) % 3] == v.a[(m + 1) % 3]) ||
								(u.a[(k + 1) % 3] == v.a[m] && u.a[k] == v.a[(m + 1) % 3])) {
								bad_edge = true;
								break;
							}
						}
						if (bad_edge) {
							break;
						}
					}
				}
				if (!bad_edge) {
					polygon.push_back({ u.a[k], u.a[(k + 1) % 3] });
				}
			}
		}
		
		triangles = gud_triangles;
		for (auto& u : polygon) {
			t.a[0] = u.first;
			t.a[1] = u.second;
			t.a[2] = i;
			triangles.push_back(t);
		}
	}
	
	gud_triangles.clear();
	for (auto& u : triangles) {
		if (u.a[0] > -1 && u.a[1] > -1 && u.a[2] > -1) {
			gud_triangles.push_back(u);
		}
	}
	return gud_triangles;
}

sorted_tri_comparator sorted_tri_pair::stc;

void t_voronoi(float* points, int N, std::vector<float>& vertices, float* border, int K,
	std::map<sorted_tri, std::vector<float>, sorted_tri_comparator>& valid_vor_verts,
	std::vector<std::vector<sorted_tri_pair>>& site_polygon,
	std::vector<std::vector<int>>& nbrdata,
	std::map<delaunay_tri_edge, std::pair<sorted_tri, sorted_tri>, delaunay_edge_comparator>& valid_edges,
	std::set<sorted_tri_pair, sorted_tri_pair_comparator>& removed_edges_set,
	int& start_ind) {
	float mmx[2];
	float mmy[2];
	get_min_max_pos(mmx, mmy, points, N, 0);
	printf("mmx: (%f, %f); mmy: (%f, %f)\n", mmx[0], mmx[1], mmy[0], mmy[1]);

	float* sites = points;
	std::vector<float> tmp;
	if (1 && border != nullptr) {
		for (int i = 0; i < N; ++i) {
			if (is_point_inside_polygon(border, K, &points[2 * i])) {
				tmp.push_back(points[2 * i]);
				tmp.push_back(points[2 * i + 1]);
			}
			else {
				//printf("not inside: [%d] (%f, %f)\n", i, points[2 * i], points[2 * i + 1]);
			}
		}
		sites = tmp.data();
		N = tmp.size() / 2;
	}
	if (sites == nullptr) {
		printf("no points inside the cull polygon\n");
		return;
	}
	std::map<delaunay_tri_edge, tri, delaunay_edge_comparator> edge_map;
	
	float a = 5 * (mmx[1] - mmx[0] > mmy[1] - mmy[0] ? mmy[1] - mmy[0] : mmx[1] - mmx[0]);
	float p[3][2] = { {mmx[0] - a, mmy[0] - a}, { mmx[0] + 2 * a, mmy[0] - a},  { mmx[0] + a / 2.0, mmy[0] + a * 1.80277564f} };
	float* s[3];

	int triangles_count = 0;
	std::vector<tri> triangles(2*N + 1);
	std::vector<tri> bad_triangles(2*N + 1);
	
	delaunay_tri_edge dte;
	
	tri t;
	t.a[0] = -3;
	t.a[1] = -2;
	t.a[2] = -1;
	triangles[triangles_count++] = t;
	dte.set(-3, -2);
	edge_map.insert({ dte, {0, -1, -4} });
	dte.set(-2, -1);
	edge_map.insert({ dte, {0, -3, -4} });
	dte.set(-1, -3);
	edge_map.insert({ dte, {0, -2, -4} });

	int bad_triangles_count = 0;

	for (int i = 0; i < N; ++i) {
		bad_triangles_count = 0;
		for (int j = 0; j < triangles_count; ++j) {
			auto& u = triangles[j];
			float q[2];
			// mid
			for (int k = 0; k < 3; ++k) {
				if (u.a[k] < 0) {
					s[k] = p[u.a[k] + 3];
				}
				else {
					s[k] = &sites[2 * u.a[k]];
				}
			}
			float c[2];
			get_tri_circumcirlce(s[0], s[1], s[2], c);
			q[0] = s[0][0] - c[0];
			q[1] = s[0][1] - c[1];
			c[0] = c[0] - sites[2 * i];
			c[1] = c[1] - sites[2 * i + 1];
			if (sqr_len(q) > sqr_len(c)) {
				for (int k = 0; k < 3; ++k) {
					dte.set(u.a[k], u.a[(k + 1) % 3]);
					//printf("test (%d, %d) v (%d)\n", dte.a, dte.b, u.a[(k + 2)%3]);
					auto res = edge_map.find(dte);
					if (res->second.a[0] == 0) {
						res->second.a[0] = 1;
						if (res->second.a[1] == u.a[(k + 2) % 3]) {
							res->second.a[1] = -4;
						}
						else {
							res->second.a[2] = -4;
						}
					}
					else {
						//printf("erase (%d, %d)\n", dte.a, dte.b);
						edge_map.erase(res);
					}
				}
				//printf("bad2: (%d, %d, %d)\n", triangles[j].a[0], triangles[j].a[1], triangles[j].a[2]);
				bad_triangles[bad_triangles_count++] = triangles[j];
				triangles[j] = triangles[triangles_count - 1];
				--triangles_count;
				--j;
			}
		}
		
		for (int j = 0; j < bad_triangles_count; ++j) {
			for (int k = 0; k < 3; ++k) {
				dte.set(bad_triangles[j].a[k], bad_triangles[j].a[(k + 1) % 3]);
				auto q = edge_map.find(dte);
				if (q != edge_map.end()) {
					auto& u = *q;
					if (u.second.a[0] == 1) {
						//printf("e: (%d, %d), add: (%d)\n", u.first.a, u.first.b, i);
						u.second.a[0] = 0;
						if (u.second.a[1] == -4) {
							u.second.a[1] = i;
						}
						else {
							u.second.a[2] = i;
						}

						t.a[0] = u.first.a;
						t.a[1] = u.first.b;
						t.a[2] = i;
						if (point_orientation(&sites[2 * u.first.a], &sites[2 * u.first.b], &sites[2 * i]) < 0) {
							t.a[0] = u.first.a;
							t.a[1] = i;
							t.a[2] = u.first.b;
						}

						dte.set(i, u.first.a);
						auto res = edge_map.find(dte);
						if (res == edge_map.end())
						{
							edge_map.insert(res, { dte, { 0, u.first.b, -4} });
						}
						else {
							res->second.a[2] = u.first.b;
						}
						dte.set(u.first.b, i);
						res = edge_map.find(dte);
						if (res == edge_map.end())
						{
							edge_map.insert(res, { dte, { 0, u.first.a, -4} });
						}
						else {
							res->second.a[2] = u.first.a;
						}
						triangles[triangles_count++] = t;
					}
				}
			}
			
		}

	}

	float min_len_threshold = 0.01;

	float f[3];
	f[0] = (rand() % 1000) / 1000.f;
	f[1] = (1 - f[0]) * ((rand() % 1000) / 1000.f);
	f[2] = 1 - f[0] - f[1];

	enum class PointSelectionOption {
		Circumcenter,
		Centeroid,
		Incenter,
		MidNinePointCenterCircumcenter,
		RandomInside
	};

	PointSelectionOption point_selection_option = PointSelectionOption::MidNinePointCenterCircumcenter;
	//std::map<sorted_tri, std::vector<float>, sorted_tri_comparator> valid_vor_verts;
	//std::map<delaunay_tri_edge, std::pair<sorted_tri, sorted_tri>, delaunay_edge_comparator> valid_edges;
	std::map<sorted_tri, sorted_tri, sorted_tri_comparator> vertex_aliases_map;
	std::vector<std::pair<std::vector<float>, std::vector<sorted_tri>>> vertex_clusters;
	std::map<sorted_tri, int, sorted_tri_comparator> in_cluster_map;
	std::set<int> has_invalid_vert_set;
	for (auto& u : edge_map) {
		if (u.first.a > -1 && u.first.b > -1) {
			float* s[2] = { &sites[2 * u.first.a], &sites[2 * u.first.b] };
			float vd; 
			float v[2];
			float c1[2];
			float c2[2];

			int c = u.second.a[1];
			int d = u.second.a[2];
			if (c > -1 && d > -1) {
				switch (point_selection_option)
				{
				case PointSelectionOption::Circumcenter: {
					get_tri_circumcirlce(s[0], s[1], &sites[2 * c], c1);
					get_tri_circumcirlce(s[0], s[1], &sites[2 * d], c2);
				} break;
				case PointSelectionOption::Centeroid: {
					c1[0] = (s[0][0] + s[1][0] + sites[2 * c]) / 3;
					c1[1] = (s[0][1] + s[1][1] + sites[2 * c + 1]) / 3;
					c2[0] = (s[0][0] + s[1][0] + sites[2 * d]) / 3;
					c2[1] = (s[0][1] + s[1][1] + sites[2 * d + 1]) / 3;
				} break;
				case PointSelectionOption::Incenter: {
					float l = sqrt(sqr_dist(s[0], s[1]));
					float l1[2] = { sqrt(sqr_dist(s[0], &sites[2 * c])), sqrt(sqr_dist(s[1], &sites[2 * c])) };
					float l2[2] = { sqrt(sqr_dist(s[0], &sites[2 * d])), sqrt(sqr_dist(s[1], &sites[2 * d])) };
					
					c1[0] = (l * sites[2 * c] + l1[0] * s[1][0] + l1[1] * s[0][0]) / (l + l1[0] + l1[1]);
					c1[1] = (l * sites[2 * c + 1] + l1[0] * s[1][1] + l1[1] * s[0][1]) / (l + l1[0] + l1[1]);

					c2[0] = (l * sites[2 * d] + l2[0] * s[1][0] + l2[1] * s[0][0]) / (l + l2[0] + l2[1]);
					c2[1] = (l * sites[2 * d + 1] + l2[0] * s[1][1] + l2[1] * s[0][1]) / (l + l2[0] + l2[1]);
					
				} break;
				case PointSelectionOption::MidNinePointCenterCircumcenter: {
					float c3[2];
					float mid_ab[2] = { (s[0][0] + s[1][0]) / 2.f, (s[0][1] + s[1][1]) / 2.f };
					float mid_a_[2] = { (s[0][0] + sites[2 * c]) / 2.f, (s[0][1] + sites[2 * c + 1]) / 2.f };
					float mid_b_[2] = { (s[1][0] + sites[2 * c]) / 2.f, (s[1][1] + sites[2 * c + 1]) / 2.f };
					get_tri_circumcirlce(mid_ab, mid_a_, mid_b_, c1);
					get_tri_circumcirlce(s[0], s[1], &sites[2*c], c3);
					c1[0] = (c1[0] + c3[0]) / 2;
					c1[1] = (c1[1] + c3[1]) / 2;

					mid_a_[0] = (s[0][0] + sites[2 * d]) / 2.f;
					mid_a_[1] = (s[0][1] + sites[2 * d + 1]) / 2.f;
					mid_b_[0] = (s[1][0] + sites[2 * d]) / 2.f;
					mid_b_[1] = (s[1][1] + sites[2 * d + 1]) / 2.f;
					get_tri_circumcirlce(mid_ab, mid_a_, mid_b_, c2);
					get_tri_circumcirlce(s[0], s[1], &sites[2 * d], c3);
					c2[0] = (c2[0] + c3[0]) / 2;
					c2[1] = (c2[1] + c3[1]) / 2;
				} break;
				case PointSelectionOption::RandomInside: {
					c1[0] = f[0] * s[0][0] + f[1] * s[1][0] + f[2] * sites[2 * c];
					c1[1] = f[0] * s[0][1] + f[1] * s[1][1] + f[2] * sites[2 * c + 1];
					c2[0] = f[0] * s[0][0] + f[1] * s[1][0] + f[2] * sites[2 * d];
					c2[1] = f[0] * s[0][1] + f[1] * s[1][1] + f[2] * sites[2 * d + 1];
				} break;
				default:
					break;
				}
			}
			else {
				for (int i = 0; i < 2; ++i)
					if (u.second.a[i + 1] < 0) {
						c = u.second.a[(i + 1) % 2 + 1];
						d = u.second.a[i + 1];
						switch (point_selection_option)
						{
						case PointSelectionOption::Circumcenter: {
							get_tri_circumcirlce(s[0], s[1], &sites[2 * c], c1);
						} break;
						case PointSelectionOption::Centeroid: {
							c1[0] = (s[0][0] + s[1][0] + sites[2 * c]) / 3;
							c1[1] = (s[0][1] + s[1][1] + sites[2 * c + 1]) / 3;
						} break;
						case PointSelectionOption::Incenter: {
							float l = sqrt(sqr_dist(s[0], s[1]));
							float l1[2] = { sqrt(sqr_dist(s[0], &sites[2 * c])), sqrt(sqr_dist(s[1], &sites[2 * c])) };
							c1[0] = (l * sites[2 * c] + l1[0] * s[1][0] + l1[1] * s[0][0]) / (l + l1[0] + l1[1]);
							c1[1] = (l * sites[2 * c + 1] + l1[0] * s[1][1] + l1[1] * s[0][1]) / (l + l1[0] + l1[1]);
						} break;
						case PointSelectionOption::MidNinePointCenterCircumcenter: {
							float mid_ab[2] = { (s[0][0] + s[1][0]) / 2.f, (s[0][1] + s[1][1]) / 2.f };
							float mid_a_[2] = { (s[0][0] + sites[2 * c]) / 2.f, (s[0][1] + sites[2 * c + 1]) / 2.f };
							float mid_b_[2] = { (s[1][0] + sites[2 * c]) / 2.f, (s[1][1] + sites[2 * c + 1]) / 2.f };
							get_tri_circumcirlce(mid_ab, mid_a_, mid_b_, c1);
							get_tri_circumcirlce(s[0], s[1], &sites[2 * c], c2);
							c1[0] = (c1[0] + c2[0]) / 2;
							c1[1] = (c1[1] + c2[1]) / 2;
							} break;
						case PointSelectionOption::RandomInside: {
							c1[0] = f[0] * s[0][0] + f[1] * s[1][0] + f[2] * sites[2 * c];
							c1[1] = f[0] * s[0][1] + f[1] * s[1][1] + f[2] * sites[2 * c + 1];
						} break;
						default:
							break;
						}

						float sign = -1;
						if (point_orientation(s[0], s[1], &sites[2 * u.second.a[(i + 1) % 2 + 1]]) < 0) {
							sign = 1;
						}

						c2[0] = sign * (s[1][1] - s[0][1]);
						c2[1] = sign * (s[0][0] - s[1][0]);
						float l = std::sqrt(sqr_len(c2));
						c2[0] *= .725f / l;
						c2[1] *= .725f / l;
						c2[0] += c1[0];
						c2[1] += c1[1];
					}
			}
			// Note/TODO: the border will be assumed to be convex in the calcs below to make calcs simpler for now
			// In the future make it work for any simple polygon
			// A possible algorithm outline 
			// 1. Get intersected parts of a 'voronoi edge' that lie inside the border
			// 2. find a shortest path from 1 site to the other
			//	. either use, in some clever way, delaunay trigangles or triangulate the polygon 
			// find the segment that intersects the shortest path - and that's our wall
			
#if 0
			bool b1 = is_point_inside_polygon(border, K, c1);
			bool b2 = is_point_inside_polygon(border, K, c2);
			if (!b1 && !b2) {
				//continue;
				float ppj[2] = { (s[0][0] + s[1][0]) / 2.f, (s[0][1] + s[1][1]) / 2.f };
				float s = -1;
				if (c2[0] == c1[0]) {
					s = (ppj[1] - c1[1]) / (c2[1] - c1[1]);
				}
				else {
					s = (ppj[0] - c1[0]) / (c2[0] - c1[0]);
				}
				if (s > 0 && s < 1 && is_point_inside_polygon(border, K, ppj)) {
					float res[2];
					float cut[2][2];
					float s1[2] = {1, 1}, s2[2];
					for (int k = 0; k < K; ++k) {
						float* pb[2] = { &border[2 * k], &border[2 * ((k + 1) % K)] };
						if (line_seg_and_seg_intesection(ppj, c1, pb[0], pb[1], res)) {
							if (ppj[0] == c1[0]) {
								s2[0] = (res[1] - ppj[1]) / (c1[1] - ppj[1]);
							}
							else {
								s2[0] = (res[0] - ppj[0]) / (c1[0] - ppj[0]);
							}
							if (s2[0] < s1[0]) {
								s1[0] = s2[0];
								cut[0][0] = res[0];
								cut[0][1] = res[1];
							}
						}

						if (line_seg_and_seg_intesection(ppj, c2, pb[0], pb[1], res)) {
							if (ppj[0] == c2[0]) {
								s2[1] = (res[1] - ppj[1]) / (c2[1] - ppj[1]);
							}
							else {
								s2[1] = (res[0] - ppj[0]) / (c2[0] - ppj[0]);
							}
							if (s2[1] < s1[1]) {
								s1[1] = s2[1];
								cut[1][0] = res[0];
								cut[1][1] = res[1];
							}
						}
					}
					c1[0] = cut[0][0];
					c1[1] = cut[0][1];
					c2[0] = cut[1][0];
					c2[1] = cut[1][1];
				}
				else {
					continue;
				}
				
			}else if (!b1 || !b2) {
				float* isd_p = c1;
				float* osd_p = c2;
				float res[2];
				float cut[2];
				float r1 = 1;
				if (b2) {
					int r = c;
					c = d;
					d = r;
					isd_p = c2;
					osd_p = c1;
				}
				for (int k = 0; k < K; ++k) {
					float* pb[2] = { &border[2 * k], &border[2 * ((k + 1) % K)] };
					if (line_seg_and_seg_intesection(c1, c2, pb[0], pb[1], res)) {
						float r2 = -1;
						if (fabs(c1[0] - c2[0]) < 0.0001) {
							r2 = (res[1] - isd_p[1]) / (osd_p[1] - isd_p[1]);
						}
						else {
							r2 = (res[0] - isd_p[0]) / (osd_p[0] - isd_p[0]);
						}
						if (r2 != -1 && r2 < r1) {
							r1 = r2;
							cut[0] = res[0];
							cut[1] = res[1];
						}
					}
				}
				c1[0] = isd_p[0];
				c1[1] = isd_p[1];
				c2[0] = cut[0];
				c2[1] = cut[1];
			}
#endif
			bool b1 = is_point_inside_polygon(border, K, c1);
			bool b2 = is_point_inside_polygon(border, K, c2);

			if (0 || b1 && b2) {
				if (0 || c > -1 && d > -1) {
					sorted_tri st1(u.first.a, u.first.b, c);
					sorted_tri st2(u.first.a, u.first.b, d);
					float mlt2 = min_len_threshold * min_len_threshold;
					if (sqr_dist(c1, c2) < mlt2) {
						auto it1 = in_cluster_map.find(st1);
						if (it1 != in_cluster_map.end()) {
							auto it2 = in_cluster_map.find(st2);
							if (it2 != in_cluster_map.end()) {
								// merge 2 clusters
								if (sqr_dist(vertex_clusters[it1->second].first.data(), vertex_clusters[it2->second].first.data()) < mlt2) {
									float M1 = vertex_clusters[it1->second].second.size();
									float M2 = vertex_clusters[it2->second].second.size();
									vertex_clusters[it1->second].first[0] =
										(vertex_clusters[it1->second].first[0] * M1 + vertex_clusters[it2->second].first[0] * M2) / (M1 + M2);
									vertex_clusters[it1->second].first[1] =
										(vertex_clusters[it1->second].first[1] * M1 + vertex_clusters[it2->second].first[1] * M2) / (M1 + M2);

									for (auto& u : vertex_clusters[it2->second].second) {
										in_cluster_map[u] = it1->second;
										vertex_clusters[it1->second].second.push_back(u);
									}
									for (auto& u : vertex_clusters.back().second) {
										in_cluster_map[u] = it2->second;
									}
									vertex_clusters[it2->second] = vertex_clusters.back();
									vertex_clusters.erase(vertex_clusters.end() - 1);
									//vertex_clusters.pop_back();
								}
								else {
									// just add as a valid edge
									valid_edges.insert({ u.first, {st1, st2} });
									valid_vor_verts.insert({ st1, {c1[0], c1[1]} });
									valid_vor_verts.insert({ st2, {c2[0], c2[1]} });
								}
							}
							else {
								if (sqr_dist(vertex_clusters[it1->second].first.data(), c2) < mlt2) {
									// add to the cluster 
									float M = vertex_clusters[it1->second].second.size();
									vertex_clusters[it1->second].first[0] =
										(vertex_clusters[it1->second].first[0] * M + c2[0]) / (M + 1);
									vertex_clusters[it1->second].first[1] =
										(vertex_clusters[it1->second].first[1] * M + c2[1]) / (M + 1);
									vertex_clusters[it1->second].second.push_back(st2);
									in_cluster_map.insert({ st2, it1->second });
								}
								else {
									// just add as a valid edge
									valid_edges.insert({ u.first, {st1, st2} });
									valid_vor_verts.insert({ st1, {c1[0], c1[1]} });
									valid_vor_verts.insert({ st2, {c2[0], c2[1]} });
								}
							}
						}
						else {
							auto it2 = in_cluster_map.find(st2);
							if (it2 != in_cluster_map.end()) {
								if (sqr_dist(vertex_clusters[it2->second].first.data(), c1) < mlt2) {
									// add to the cluster 
									float M = vertex_clusters[it2->second].second.size();
									vertex_clusters[it2->second].first[0] =
										(vertex_clusters[it2->second].first[0] * M + c1[0]) / (M + 1);
									vertex_clusters[it2->second].first[1] =
										(vertex_clusters[it2->second].first[1] * M + c1[1]) / (M + 1);
									vertex_clusters[it2->second].second.push_back(st1);
									in_cluster_map.insert({ st1, it2->second });
								}
								else {
									// just add as a valid edge
									valid_edges.insert({ u.first, {st1, st2} });
									valid_vor_verts.insert({ st1, {c1[0], c1[1]} });
									valid_vor_verts.insert({ st2, {c2[0], c2[1]} });
								}
							}
							else {
								int ind = -1;
								float mind2 = 1e9;
								float mid[2] = { (c1[0] + c2[0]) / 2, (c1[1] + c2[1]) / 2 };
								for (int i = 0; i < vertex_clusters.size(); ++i) {
									auto& u = vertex_clusters[i];
									if (sqr_dist(mid, u.first.data()) < mind2) {
										ind = i;
										mind2 = sqr_dist(mid, u.first.data());
									}
								}
								if (ind != -1 && mind2 < mlt2) {
									float M = vertex_clusters[ind].second.size();
									vertex_clusters[ind].first[0] = (vertex_clusters[ind].first[0] * M + c1[0] + c2[0]) / (M + 2);
									vertex_clusters[ind].first[1] = (vertex_clusters[ind].first[1] * M + c1[1] + c2[1]) / (M + 2);
									vertex_clusters[ind].second.push_back(st1);
									vertex_clusters[ind].second.push_back(st2);
									in_cluster_map.insert({ st1, ind });
									in_cluster_map.insert({ st2, ind });
								}
								else {
									// create a new cluster
									vertex_clusters.push_back({ {(c1[0] + c2[0]) / 2, (c1[1] + c2[1]) / 2}, {st1, st2} });
									in_cluster_map.insert({ st1, vertex_clusters.size() - 1 });
									in_cluster_map.insert({ st2, vertex_clusters.size() - 1 });
								}
							}
						}
					}
					else {
						valid_edges.insert({ u.first, {st1, st2} });
						valid_vor_verts.insert({ st1, {c1[0], c1[1]} });
						valid_vor_verts.insert({ st2, {c2[0], c2[1]} });
					}
				}
				else {
					has_invalid_vert_set.insert(u.first.a);
					has_invalid_vert_set.insert(u.first.b);
				}
			}
			else {
				has_invalid_vert_set.insert(u.first.a);
				has_invalid_vert_set.insert(u.first.b);
			}
		}
	}
	
	printf("number of vertex clusters: %lld\n", vertex_clusters.size());
	for (int i = 0; i < vertex_clusters.size(); ++i) {
		auto& u = vertex_clusters[i];
		int ind = 0;
		valid_vor_verts[u.second[ind]] = u.first;
		for (int j = 0; j < u.second.size(); ++j) {
			if (j == ind) continue;
			vertex_aliases_map[u.second[j]] = u.second[ind];
		}
	}
	site_polygon.resize(N);
	//std::vector<std::vector<sorted_tri_pair>> site_polygon(N);

	std::map<sorted_tri, std::pair<std::vector<sorted_tri>, std::vector<std::vector<float>>>, sorted_tri_comparator> pmap;
	//int start_ind;
	nbrdata.resize(N);
	for (auto& u : valid_edges) {	
		auto p1 = u.second.first;
		auto p2 = u.second.second;

		auto it2 = vertex_aliases_map.find(p1);
		if (it2 != vertex_aliases_map.end()) {
			p1 = it2->second;
		}
		it2 = vertex_aliases_map.find(p2);
		if (it2 != vertex_aliases_map.end()) {
			p2 = it2->second;
		}

		u.second.first = p1;
		u.second.second = p2;

		auto it = pmap.find(p1);
		if (it == pmap.end()) {
			pmap.insert(it, { p1, {{p2}, {}} });
		}
		else {
			it->second.first.push_back(p2);
		}
		it = pmap.find(p2);
		if (it == pmap.end()) {
			pmap.insert(it, { p2, {{p1}, {}} });
		}
		else {
			it->second.first.push_back(p1);
		}
		
		sorted_tri_pair stp(p1, p2);
		if (has_invalid_vert_set.find(u.first.a) == has_invalid_vert_set.end()) {
			site_polygon[u.first.a].push_back(stp);
		}
		if (has_invalid_vert_set.find(u.first.b) == has_invalid_vert_set.end()) {
			site_polygon[u.first.b].push_back(stp);
		}

		if (has_invalid_vert_set.find(u.first.a) == has_invalid_vert_set.end() && has_invalid_vert_set.find(u.first.b) == has_invalid_vert_set.end()) {
			start_ind = u.first.a;
			nbrdata[u.first.a].push_back(u.first.b);
			nbrdata[u.first.b].push_back(u.first.a);
			dte.set(u.first.a, u.first.b);
		}
		
	}
	
	while (true) {
		bool found = false;
		auto it = pmap.begin();
		for (; it != pmap.end(); ++it) {
			if (it->second.first.size() == 1) {
				found = true;
				break;
			}
		}
		if (found) {
			sorted_tri cand = it->first;
			sorted_tri other = it->second.first[0];
			pmap.erase(it);
			auto& nlist = pmap[other].first;
			auto it2 = nlist.begin();
			for (; it2 != nlist.end(); ++it2) {
				if (it2->a[0] == cand.a[0] && it2->a[1] == cand.a[1] && it2->a[2] == cand.a[2]) {
					break;
				}
			}
			nlist.erase(it2);
		}
		else {
			break;
		}
	}

	//std::set<sorted_tri_pair, sorted_tri_pair_comparator> removed_n_edges;
	std::vector<bool> is_visited(N, false);
	std::vector<int> cells;
	
	cells.push_back(start_ind);
	is_visited[start_ind] = true;
	while (1 && !cells.empty()) {
		int ind = cells.size() - 1;
		if (rand() % 2 == 1) {
			ind = rand() % cells.size();
		}
		int c = cells[ind];
		cells[ind] = cells.back();
		cells.erase(cells.end() - 1);
		for (int i = 0; i < nbrdata[c].size(); ++i) {
			if (!is_visited[nbrdata[c][i]]) {
				cells.push_back(c);
				is_visited[nbrdata[c][i]] = true;
				cells.push_back(nbrdata[c][i]);

				dte.set(c, nbrdata[c][i]);
				auto it = valid_edges.find(dte);
				sorted_tri& st1 = it->second.first;
				sorted_tri& st2 = it->second.second;
				

				/*auto it2 = vertex_aliases_map.find(st1);
				if (it2 != vertex_aliases_map.end()) {
					st1 = it2->second;
				}
				it2 = vertex_aliases_map.find(st2);
				if (it2 != vertex_aliases_map.end()) {
					st2 = it2->second;
				}*/

				auto* nlist = &pmap[st1].first;
				auto it3 = nlist->begin();
				for(;it3 != nlist->end(); ++it3) {
					if (it3->a[0] == st2.a[0] && it3->a[1] == st2.a[1] && it3->a[2] == st2.a[2]) {
						break;
					}
				}
				nlist->erase(it3);

				nlist = &pmap[st2].first;
				it3 = nlist->begin();
				for (; it3 != nlist->end(); ++it3) {
					if (it3->a[0] == st1.a[0] && it3->a[1] == st1.a[1] && it3->a[2] == st1.a[2]) {
						break;
					}
				}
				nlist->erase(it3);

				//valid_edges.erase(dte);

				sorted_tri_pair stp(st1, st2);
				removed_edges_set.insert(stp);
				break;
			}
		}
	}
	//start_ind = 1303;
	printf("start ind: %d\n", start_ind);

	float hfw = 0.0065;

	for (auto& u : pmap) {
		float* center = valid_vor_verts[u.first].data();
		int M = u.second.first.size();
		u.second.second.resize(M);
		std::vector<float> indices(M);
		std::vector<float> pos;
		
		for (int i = 0; i < M; ++i) {
			indices[i] = i;
			float* v = valid_vor_verts[u.second.first[i]].data();
			pos.push_back(v[0] - center[0]);
			pos.push_back(v[1] - center[1]);
		}
		
		// perform vec sort
		if (M > 3) {
			float o[] = { 0, 0 };

			// split points that lie on the left and right of the point[0]
			int left_pos = 1;
			int right_pos = M - 1;

			while (left_pos <= right_pos) {
				bool b1 = point_orientation(o, &pos[0], &pos[2 * left_pos]) > 0;
				bool b2 = point_orientation(o, &pos[0], &pos[2 * right_pos]) <= 0;
				if (b1 && b2) {
					//swap
					float x = pos[2 * left_pos];
					float y = pos[2 * left_pos + 1];
					pos[2 * left_pos] = pos[2 * right_pos];
					pos[2 * left_pos + 1] = pos[2 * right_pos + 1];
					pos[2 * right_pos] = x;
					pos[2 * right_pos + 1] = y;

					int t = indices[left_pos];
					indices[left_pos] = indices[right_pos];
					indices[right_pos] = t;

					++left_pos;
					--right_pos;

				}
				else {
					if (!b1) {
						++left_pos;
					}
					else {
						--right_pos;
					}
				}
			}
			//printf("left size: %d right size: %d\n", (left_pos - 1), (N - 1 - right_pos));
			if (left_pos - 1 > 1) {
				int next = 1;
				while (next != left_pos) {
					int ind = next;
					for (int i = next + 1; i < left_pos; ++i) {
						if (point_orientation(o, &pos[2 * ind], &pos[2 * i]) > 0) {
							ind = i;
						}
					}
					{
						float x = pos[2 * next];
						float y = pos[2 * next + 1];
						pos[2 * next] = pos[2 * ind];
						pos[2 * next + 1] = pos[2 * ind + 1];
						pos[2 * ind] = x;
						pos[2 * ind + 1] = y;

						int t = indices[next];
						indices[next] = indices[ind];
						indices[ind] = t;
					}
					++next;
				}
			}
			if (M - 1 - right_pos > 1) {
				int next = M - 1;
				while (next != right_pos) {
					int ind = next;
					for (int i = next - 1; i > right_pos; --i) {
						if (point_orientation(o, &pos[2 * ind], &pos[2 * i]) < 0) {
							ind = i;
						}
					}
					{
						float x = pos[2 * next];
						float y = pos[2 * next + 1];
						pos[2 * next] = pos[2 * ind];
						pos[2 * next + 1] = pos[2 * ind + 1];
						pos[2 * ind] = x;
						pos[2 * ind + 1] = y;

						int t = indices[next];
						indices[next] = indices[ind];
						indices[ind] = t;
					}
					--next;
				}
			}

			std::vector<sorted_tri> tmp(M);
			for (int i = 0; i < M; ++i) {
				tmp[i] = u.second.first[indices[i]];
			}
			u.second.first = tmp;
		}
		if (M > 2) 
		{
			auto& nlist = u.second.first;
			double sa = 0;
			for (int i = 0; i < M - 1; ++i) {
				float* p0 = valid_vor_verts[nlist[i]].data();
				float* p1 = valid_vor_verts[nlist[i + 1]].data();
				sa += (p1[0] - p0[0]) * (p1[1] + p0[1]);
			}
			float* p0 = valid_vor_verts[nlist[0]].data();
			float* p1 = valid_vor_verts[nlist.back()].data();
			sa += (p0[0] - p1[0]) * (p1[1] + p0[1]);

			if (sa > 0) {
				for (auto li = nlist.begin(), ri = nlist.end() - 1; li <= ri; ++li, --ri) {
					auto tp = *li;
					*li = *ri;
					*ri = tp;
				}

				//printf("before CW\n");
			}
			else {
				//printf("before CCW\n");
			}
		}
	}

	
	auto get_bisector_vec_with_mag = [](float* p1, float* p2, float* p3, float hfw, float* out) {
		float a = sqrt(sqr_dist(p2, p3));
		float b = sqrt(sqr_dist(p1, p3));
		float c = sqrt(sqr_dist(p1, p2));
		float av[2] = { p2[0] - p3[0], p2[1] - p3[1] };
		float bv[2] = { p1[0] - p3[0], p1[1] - p3[1] };
		float d = b / (b + c);
		float p[2] = { av[0] * d + p3[0], av[1] * d + p3[1] };
		
		float h = dot_prd(av, bv) / (b + c);
		h = sqrt(d * d * a * a - h * h);
		float mag = hfw / h;
		out[0] = mag * (p[0] - p1[0]);
		out[1] = mag * (p[1] - p1[1]);

		/*float av[2] = { p2[0] - p1[0], p2[1] - p1[1] };
		float bv[2] = { p3[0] - p1[0], p3[1] - p1[1] };
		float l1 = sqrt(sqr_len(av));
		float l2 = sqrt(sqr_len(bv));
		out[0] = av[0] * l2 + bv[0] * l1;
		out[1] = av[1] * l2 + bv[1] * l1;
		l1 = sqrt(sqr_len(out));
		out[0] *= hfw * 2.2 / l1;
		out[1] *= hfw * 2.2 / l1;*/
	};
		 
	std::set<sorted_tri_pair, sorted_tri_pair_comparator> done_set;
	sorted_tri_pair stp;
#if 0
	auto get_vs = [&get_bisector_vec_with_mag](float* p1, float* p2, float* p1prev, float* p1next, float hfw, float* u1, float* u2) {
		if (p1prev == nullptr) {
			// only 1 edge, N = 1
			u1[0] = p1[1] - p2[1];
			u1[1] = p2[0] - p1[0];
			float d = sqrt(sqr_len(u1));
			u1[0] *= hfw / d;
			u1[1] *= hfw / d;
			u2[0] = -u1[0];
			u2[1] = -u1[1];
			return;
		}
		if (p1next == nullptr) {
			// only 2 edge, N = 2
			float q1[2] = { p1[0] - p1prev[0], p1[1] - p1prev[1] };
			float q2[2] = { p2[0] - p1[0], p2[1] - p1[1] };
			float d = dot_prd(q1, q2);
			float th = 0.99;
			if (d * d > th * th * sqr_len(q1) * sqr_len(q2)) {
				u1[0] = p1[1] - p2[1];
				u1[1] = p2[0] - p1[0];
				float d = sqrt(sqr_len(u1));
				u1[0] *= hfw / d;
				u1[1] *= hfw / d;
				u2[0] = -u1[0];
				u2[1] = -u1[1];
				return;
			}
			get_bisector_vec_with_mag(p1, p1prev, p2, hfw, u1);
			u2[0] = -u1[0];
			u2[1] = -u1[1];
			return;
		}

		get_bisector_vec_with_mag(p1, p1prev, p2, hfw, u1);
		get_bisector_vec_with_mag(p1, p1next, p2, hfw, u2);

		float s1 = point_orientation(p1, p2, p1prev);
		float s2 = point_orientation(p1, p2, p1next);

		// if both lie on the same side
		if ((s1 > 0 && s2 > 0) || (s1 < 0 && s2 < 0)) {
			float s3 = point_orientation(p1, p1prev, p1next);
			float s4 = point_orientation(p1, p1prev, p2);
			if ((s3 > 0 && s4 < 0) || (s3 < 0 && s4 > 0)) {
				u2[0] = -u2[0];
				u2[1] = -u2[1];
			}
			else {
				u1[0] = -u1[0];
				u1[1] = -u1[1];
			}
		}
	};
	for (auto& u : pmap) {
		float* oo = valid_vor_verts[u.first].data();
		float* pp = oo;
		auto& nlist = u.second.first;
		int M1 = nlist.size();
		for (int i = 0; i < M1; ++i) {
			stp.set(u.first, nlist[i]);
			auto it = done_set.find(stp);
			if (it != done_set.end()) {
				continue;
			}
			done_set.insert(it, stp);

			float* tt = valid_vor_verts[nlist[i]].data();
			float dx = pp[0] - tt[0];
			float dy = pp[1] - tt[1];
			{
				float u1[2] = { 0, 0 };
				float u2[2] = { 0, 0 };

				if (M1 == 1) {
					get_vs(pp, tt, nullptr, nullptr, hfw, u1, u2);
				}
				else if (M1 == 2) {
					get_vs(pp, tt, valid_vor_verts[nlist[(i - 1 + M1) % M1]].data(), nullptr, hfw, u1, u2);
				}
				else {
					get_vs(pp, tt, valid_vor_verts[nlist[(i - 1 + M1) % M1]].data(), valid_vor_verts[nlist[(i + 1) % M1]].data(), hfw, u1, u2);
				}

				auto& mlist = pmap[nlist[i]].first;
				int M2 = mlist.size();
				for (int j = 0; j < M2; ++j) {
					if (u.first.a[0] == mlist[j].a[0] && u.first.a[1] == mlist[j].a[1] && u.first.a[2] == mlist[j].a[2]) {
						float v1[2] = { 0, 0 };
						float v2[2] = { 0, 0 };

						if (M2 == 1) {
							get_vs(tt, pp, nullptr, nullptr, hfw, v1, v2);
						}
						else if (M2 == 2) {
							get_vs(tt, pp, valid_vor_verts[mlist[(j - 1 + M2) % M2]].data(), nullptr, hfw, v1, v2);
						}
						else {
							get_vs(tt, pp, valid_vor_verts[mlist[(j - 1 + M2) % M2]].data(), valid_vor_verts[mlist[(j + 1) % M2]].data(), hfw, v1, v2);
						}

						float l1[2] = { pp[0] + u1[0], pp[1] + u1[1] };
						float l2[2] = { tt[0] + v1[0], tt[1] + v1[1] };
						float l3[2] = { pp[0] + u2[0], pp[1] + u2[1] };
						float l4[2] = { tt[0] + v2[0], tt[1] + v2[1] };
						float* w1 = v1, * w2 = v2;
						float res[2];
						if (line_seg_and_seg_intesection(l1, l2, l3, l4, res)) {
							w1 = v2;
							w2 = v1;
						}
						if (0) {
							vertices.push_back(pp[0] + u1[0]);
							vertices.push_back(pp[1] + u1[1]);
							vertices.push_back(tt[0] + w1[0]);
							vertices.push_back(tt[1] + w1[1]);
							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);

							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);
							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
							vertices.push_back(pp[0] + u1[0]);
							vertices.push_back(pp[1] + u1[1]);

							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);
							vertices.push_back(tt[0] + w2[0]);
							vertices.push_back(tt[1] + w2[1]);
							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);

							vertices.push_back(tt[0] + w2[0]);
							vertices.push_back(tt[1] + w2[1]);
							vertices.push_back(pp[0] + u2[0]);
							vertices.push_back(pp[1] + u2[1]);
							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
						}
						else if (1) {
							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);

							vertices.push_back(pp[0] + u1[0]);
							vertices.push_back(pp[1] + u1[1]);
							vertices.push_back(tt[0] + w1[0]);
							vertices.push_back(tt[1] + w1[1]);

							vertices.push_back(tt[0] + w1[0]);
							vertices.push_back(tt[1] + w1[1]);
							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);

							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);
							vertices.push_back(pp[0] + u1[0]);
							vertices.push_back(pp[1] + u1[1]);

							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);
							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);

							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
							vertices.push_back(pp[0] + u1[0]);
							vertices.push_back(pp[1] + u1[1]);

							vertices.push_back(pp[0] + u1[0]);
							vertices.push_back(pp[1] + u1[1]);
							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);


							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);
							vertices.push_back(tt[0] + w2[0]);
							vertices.push_back(tt[1] + w2[1]);

							vertices.push_back(tt[0] + w2[0]);
							vertices.push_back(tt[1] + w2[1]);
							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);

							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);

							vertices.push_back(tt[0] + w2[0]);
							vertices.push_back(tt[1] + w2[1]);
							vertices.push_back(pp[0] + u2[0]);
							vertices.push_back(pp[1] + u2[1]);

							vertices.push_back(pp[0] + u2[0]);
							vertices.push_back(pp[1] + u2[1]);
							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);

							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
							vertices.push_back(tt[0] + w2[0]);
							vertices.push_back(tt[1] + w2[1]);
						}
						else {
							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);

							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
							vertices.push_back(pp[0] + u1[0]);
							vertices.push_back(pp[1] + u1[1]);

							vertices.push_back(pp[0]);
							vertices.push_back(pp[1]);
							vertices.push_back(pp[0] + u2[0]);
							vertices.push_back(pp[1] + u2[1]);

							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);
							vertices.push_back(tt[0] + w1[0]);
							vertices.push_back(tt[1] + w1[1]);

							vertices.push_back(tt[0]);
							vertices.push_back(tt[1]);
							vertices.push_back(tt[0] + w2[0]);
							vertices.push_back(tt[1] + w2[1]);
						}
						break;
					}
				}
			}
		}
	}
#else
	auto get_vs2 = [&get_bisector_vec_with_mag](float* p1, float* p2, float* p1prev, float* p1next, float hfw, float* u1, float* u2) {
		if (p1prev == nullptr) {
			// only 1 edge, N = 1
			u1[0] = p2[1] - p1[1];
			u1[1] = p1[0] - p2[0];
			float d = sqrt(sqr_len(u1));
			u1[0] *= hfw / d;
			u1[1] *= hfw / d;
			u2[0] = -u1[0];
			u2[1] = -u1[1];
			return;
		}
		if (p1next == nullptr) {
			// only 2 edge, N = 2
			float q1[2] = { p1[0] - p1prev[0], p1[1] - p1prev[1] };
			float q2[2] = { p2[0] - p1[0], p2[1] - p1[1] };
			float d = dot_prd(q1, q2);
			float th = 0.99;
			if (d * d > th * th * sqr_len(q1) * sqr_len(q2)) {
				u1[0] = p1[1] - p2[1];
				u1[1] = p2[0] - p1[0];
				float d = sqrt(sqr_len(u1));
				u1[0] *= hfw / d;
				u1[1] *= hfw / d;
				u2[0] = -u1[0];
				u2[1] = -u1[1];
			}
			get_bisector_vec_with_mag(p1, p1prev, p2, hfw, u1);
			q1[0] = p1[0] + u1[0];
			q1[1] = p1[1] + u1[1];
			if (point_orientation(p1, p2, q1) < 0) {
				u1[0] *= -1;
				u1[1] *= -1;
			}
			u2[0] = -u1[0];
			u2[1] = -u1[1];
			return;
		}

		get_bisector_vec_with_mag(p1, p1prev, p2, hfw, u1);
		get_bisector_vec_with_mag(p1, p1next, p2, hfw, u2);

		float q[2] = { u1[0] + p1[0], u1[1] + p1[1] };
		if (point_orientation(p1, p2, q) < 0) {
			u1[0] *= -1;
			u1[1] *= -1;
		}
		q[0] = u2[0] + p1[0];
		q[1] = u2[1] + p1[1];
		if (point_orientation(p1, p2, q) > 0) {
			u2[0] *= -1;
			u2[1] *= -1;
		}
	};
	
	bool line_mode = false;

	for (auto& u : pmap) {
		float* pp = valid_vor_verts[u.first].data();
		auto& nlist = u.second.first;
		int M1 = nlist.size();

		for (int i = 0; i < M1; ++i) {
			stp.set(u.first, nlist[i]);
			auto it = done_set.find(stp);
			if (it != done_set.end()) {
				continue;
			}
			done_set.insert(it, stp);

			float* tt = valid_vor_verts[nlist[i]].data();
			
			float u1[2] = { 0, 0 };
			float u2[2] = { 0, 0 };
			
			if (M1 == 1) {
				get_vs2(pp, tt, nullptr, nullptr, hfw, u1, u2);
			}
			else if (M1 == 2) {
				get_vs2(pp, tt, valid_vor_verts[nlist[(i - 1 + M1) % M1]].data(), nullptr, hfw, u1, u2);
			}
			else {
				get_vs2(pp, tt, valid_vor_verts[nlist[(i - 1 + M1) % M1]].data(), valid_vor_verts[nlist[(i + 1) % M1]].data(), hfw, u1, u2);
			}
			
			auto& v = pmap[nlist[i]];
			auto& mlist = v.first;
			int M2 = mlist.size();
			for (int j = 0; j < M2; ++j) {
				if (u.first.a[0] == mlist[j].a[0] && u.first.a[1] == mlist[j].a[1] && u.first.a[2] == mlist[j].a[2]) {
					float v1[2] = { 0, 0 };
					float v2[2] = { 0, 0 };

					if (M2 == 1) {
						get_vs2(tt, pp, nullptr, nullptr, hfw, v1, v2);
					}
					else if (M2 == 2) {
						get_vs2(tt, pp, valid_vor_verts[mlist[(j - 1 + M2) % M2]].data(), nullptr, hfw, v1, v2);
					}
					else {
						get_vs2(tt, pp, valid_vor_verts[mlist[(j - 1 + M2) % M2]].data(), valid_vor_verts[mlist[(j + 1) % M2]].data(), hfw, v1, v2);
					}

					//float w1 = 0.5, w2 = 1 - w1;
					float lu1 = sqrt(sqr_len(u1));
					float lu2 = sqrt(sqr_len(u2));
					float lv1 = sqrt(sqr_len(v1));
					float lv2 = sqrt(sqr_len(v2));
					float w1 = lv2 / (lu1 + lv2);
					float w2 = 1 - w1;
					float w3 = lv1 / (lu2 + lv1);
					float w4 = 1 - w3;
					
					float l = hfw / sqrt(sqr_dist(pp, tt));
					float pu1[2] = { (tt[1] - pp[1]) * l, (pp[0] - tt[0]) * l };
					float pu2[2] = { -pu1[0], -pu1[1] };

					float mid[2] = { (pp[0] + tt[0]) / 2 , (pp[1] + tt[1]) / 2 };
					u.second.second[i].push_back(u1[0] + pp[0]);
					u.second.second[i].push_back(u1[1] + pp[1]);
					//u.second.second[i].push_back(u1[0] * w1 + v2[0] * w2 + mid[0]);
					//u.second.second[i].push_back(u1[1] * w1 + v2[1] * w2 + mid[1]);
					u.second.second[i].push_back(pu1[0] + mid[0]);
					u.second.second[i].push_back(pu1[1] + mid[1]);

					u.second.second[i].push_back(u2[0] + pp[0]);
					u.second.second[i].push_back(u2[1] + pp[1]);
					//u.second.second[i].push_back(u2[0] * w3 + v1[0] * w4 + mid[0]);
					//u.second.second[i].push_back(u2[1] * w3 + v1[1] * w4 + mid[1]);
					u.second.second[i].push_back(pu2[0] + mid[0]);
					u.second.second[i].push_back(pu2[1] + mid[1]);

					v.second[j].push_back(v1[0] + tt[0]);
					v.second[j].push_back(v1[1] + tt[1]);
					//v.second[j].push_back(v1[0] * w4 + u2[0] * w3 + mid[0]);
					//v.second[j].push_back(v1[1] * w4 + u2[1] * w3 + mid[1]);
					v.second[j].push_back(pu2[0] + mid[0]);
					v.second[j].push_back(pu2[1] + mid[1]);

					v.second[j].push_back(v2[0] + tt[0]);
					v.second[j].push_back(v2[1] + tt[1]);
					//v.second[j].push_back(v2[0] * w2 + u1[0] * w1 + mid[0]);
					//v.second[j].push_back(v2[1] * w2 + u1[1] * w1 + mid[1]);
					v.second[j].push_back(pu1[0] + mid[0]);
					v.second[j].push_back(pu1[1] + mid[1]);

					if (line_mode) {
						vertices.push_back(pp[0]);
						vertices.push_back(pp[1]);
						vertices.push_back(tt[0]);
						vertices.push_back(tt[1]);

						vertices.push_back(pp[0]);
						vertices.push_back(pp[1]);
						vertices.push_back(pp[0] + u1[0]);
						vertices.push_back(pp[1] + u1[1]);

						vertices.push_back(pp[0]);
						vertices.push_back(pp[1]);
						vertices.push_back(pp[0] + u2[0]);
						vertices.push_back(pp[1] + u2[1]);

						vertices.push_back(tt[0]);
						vertices.push_back(tt[1]);
						vertices.push_back(tt[0] + v1[0]);
						vertices.push_back(tt[1] + v1[1]);

						vertices.push_back(tt[0]);
						vertices.push_back(tt[1]);
						vertices.push_back(tt[0] + v2[0]);
						vertices.push_back(tt[1] + v2[1]);
					}
					else if (0) {
						vertices.push_back(pp[0] + u1[0]);
						vertices.push_back(pp[1] + u1[1]);
						vertices.push_back(tt[0] + v2[0]);
						vertices.push_back(tt[1] + v2[1]);
						vertices.push_back(tt[0]);
						vertices.push_back(tt[1]);

						vertices.push_back(tt[0]);
						vertices.push_back(tt[1]);
						vertices.push_back(pp[0]);
						vertices.push_back(pp[1]);
						vertices.push_back(pp[0] + u1[0]);
						vertices.push_back(pp[1] + u1[1]);

						vertices.push_back(tt[0]);
						vertices.push_back(tt[1]);
						vertices.push_back(tt[0] + v1[0]);
						vertices.push_back(tt[1] + v1[1]);
						vertices.push_back(pp[0]);
						vertices.push_back(pp[1]);

						vertices.push_back(tt[0] + v1[0]);
						vertices.push_back(tt[1] + v1[1]);
						vertices.push_back(pp[0] + u2[0]);
						vertices.push_back(pp[1] + u2[1]);
						vertices.push_back(pp[0]);
						vertices.push_back(pp[1]);
					}
					break;
				}
			}
		}
	}

	// L = (a, b, c) where ax + by + c = 0
	auto get_line_and_bezier_intersection_t = [](float* L, float* p0, float* p1, float* p2) {
		float Q[3] = { dot_prd(L, p0), dot_prd(L, p1), dot_prd(L, p2) };
		float a = Q[0] - 2 * Q[1] + Q[2];
		float b = -2 * Q[0] + 2 * Q[1];
		float c = Q[0] + L[2];
		float Dt = b * b - 4 * a * c;
		float t = (-b - sqrt(Dt)) / (2 * a);
		if (t < 0 || t > 1) {
			t = (-b + sqrt(Dt)) / (2 * a);
		}
		return t;
	};

	for (auto& u : pmap) {
		float* pp = valid_vor_verts[u.first].data();
		auto& nlist = u.second.first;
		int M = nlist.size();
		if (M == 1) {
			float* A = valid_vor_verts[u.first].data();
			float* B = valid_vor_verts[u.second.first[0]].data();
			float AB[2] = { B[0] - A[0], B[1] - A[1] };
			float* p0 = &u.second.second[0][0];
			float* p3 = &u.second.second[0][4];
			
			if (line_mode) {
				vertices.push_back(p0[0]);
				vertices.push_back(p0[1]);
				vertices.push_back(u.second.second[0][2]);
				vertices.push_back(u.second.second[0][3]);

				vertices.push_back(p3[0]);
				vertices.push_back(p3[1]);
				vertices.push_back(u.second.second[0][6]);
				vertices.push_back(u.second.second[0][7]);
			}
			else {
				vertices.push_back(A[0]);
				vertices.push_back(A[1]);
				vertices.push_back(p0[0]);
				vertices.push_back(p0[1]);
				vertices.push_back(u.second.second[0][2]);
				vertices.push_back(u.second.second[0][3]);

				vertices.push_back(A[0]);
				vertices.push_back(A[1]);
				vertices.push_back(u.second.second[0][2]);
				vertices.push_back(u.second.second[0][3]);
				vertices.push_back((A[0] + B[0]) / 2);
				vertices.push_back((A[1] + B[1]) / 2);
				

				vertices.push_back(A[0]);
				vertices.push_back(A[1]);
				vertices.push_back((A[0] + B[0]) / 2);
				vertices.push_back((A[1] + B[1]) / 2);
				vertices.push_back(u.second.second[0][6]);
				vertices.push_back(u.second.second[0][7]);

				vertices.push_back(A[0]);
				vertices.push_back(A[1]);
				vertices.push_back(u.second.second[0][6]);
				vertices.push_back(u.second.second[0][7]);
				vertices.push_back(p3[0]);
				vertices.push_back(p3[1]);

			}
			

			float v[2] = { (p0[0] - p3[0]) / 2, (p0[1] - p3[1]) / 2 };
			float y = -v[0];
			v[0] = v[1];
			v[1] = y;
			if (dot_prd(AB, v) > 0) {
				v[0] *= -1;
				v[1] *= -1;
			}
			AB[0] = v[0] + A[0];
			AB[1] = v[1] + A[1];
			
			float p1[2] = { p0[0] + v[0], p0[1] + v[1] };
			float p2[2] = { p3[0] + v[0], p3[1] + v[1] };
			if(0)
			{
				vertices.push_back(p0[0]);
				vertices.push_back(p0[1]);
				vertices.push_back(p1[0]);
				vertices.push_back(p1[1]);

				vertices.push_back(p1[0]);
				vertices.push_back(p1[1]);
				vertices.push_back(p2[0]);
				vertices.push_back(p2[1]);

				vertices.push_back(p2[0]);
				vertices.push_back(p2[1]);
				vertices.push_back(p3[0]);
				vertices.push_back(p3[1]);
			}
			
			if (line_mode) {
				float step = 0.1;
				for (float t = 0; t <= 1.f - step + 0.001f; t += step) {
					float t1 = 1 - t;
					float t2 = 3 * t * t1 * t1;
					float t3 = 3 * t * t * t1;
					float t4 = t * t * t;
					t1 = t1 * t1 * t1;

					float x = t1 * p0[0] + t2 * p1[0] + t3 * p2[0] + t4 * p3[0];
					float y = t1 * p0[1] + t2 * p1[1] + t3 * p2[1] + t4 * p3[1];

					vertices.push_back(x);
					vertices.push_back(y);

					t += step;
					t1 = 1 - t;
					t2 = 3 * t * t1 * t1;
					t3 = 3 * t * t * t1;
					t4 = t * t * t;
					t1 = t1 * t1 * t1;
					t -= step;

					x = t1 * p0[0] + t2 * p1[0] + t3 * p2[0] + t4 * p3[0];
					y = t1 * p0[1] + t2 * p1[1] + t3 * p2[1] + t4 * p3[1];


					vertices.push_back(x);
					vertices.push_back(y);
				}
			}
			else {
				float step = 0.1;
				for (float t = 0; t <= 1.f - step + 0.001f; t += step) {
					float t1 = 1 - t;
					float t2 = 3 * t * t1 * t1;
					float t3 = 3 * t * t * t1;
					float t4 = t * t * t;
					t1 = t1 * t1 * t1;

					float x = t1 * p0[0] + t2 * p1[0] + t3 * p2[0] + t4 * p3[0];
					float y = t1 * p0[1] + t2 * p1[1] + t3 * p2[1] + t4 * p3[1];

					vertices.push_back(A[0]);
					vertices.push_back(A[1]);
					vertices.push_back(x);
					vertices.push_back(y);

					t += step;
					t1 = 1 - t;
					t2 = 3 * t * t1 * t1;
					t3 = 3 * t * t * t1;
					t4 = t * t * t;
					t1 = t1 * t1 * t1;
					t -= step;

					x = t1 * p0[0] + t2 * p1[0] + t3 * p2[0] + t4 * p3[0];
					y = t1 * p0[1] + t2 * p1[1] + t3 * p2[1] + t4 * p3[1];


					vertices.push_back(x);
					vertices.push_back(y);
				}
			}
		}
		else {
			for (int i = 0; i < M; ++i) {
				
				int prev_ind = (i - 1 + M) % M;
				int next_ind = (i + 1) % M;
				float* A = valid_vor_verts[u.first].data();
				float* B = valid_vor_verts[u.second.first[prev_ind]].data();
				float* C = valid_vor_verts[u.second.first[next_ind]].data();
				float* D = valid_vor_verts[u.second.first[i]].data();


				float* p0_r = &u.second.second[prev_ind][6];
				float* p1_r = &u.second.second[i][0];
				float* p2_r = &u.second.second[i][2];

				if (0) {
					vertices.push_back(p0_r[0]);
					vertices.push_back(p0_r[1]);
					vertices.push_back(p1_r[0]);
					vertices.push_back(p1_r[1]);

					vertices.push_back(p1_r[0]);
					vertices.push_back(p1_r[1]);
					vertices.push_back(p2_r[0]);
					vertices.push_back(p2_r[1]);
				}

				if (line_mode)
				{

					float step = 0.1;

					float w[3] = { 1, 1, 1 };
					float sum_w = 0;
					
					for (float t = 0; t <= 1.f - step + 0.001f; t += step) {
						float omt = 1 - t;
						sum_w = omt* omt* w[0] + 2 * omt * t * w[1] + t * t * w[2];
						float x = (omt * omt * p0_r[0] * w[0] + 2 * omt * t * p1_r[0] * w[1] + t * t * p2_r[0] * w[2]) / sum_w;
						float y = (omt * omt * p0_r[1] * w[0] + 2 * omt * t * p1_r[1] * w[1] + t * t * p2_r[1] * w[2]) / sum_w;
						vertices.push_back(x);
						vertices.push_back(y);

						t += step;
						omt = 1 - t;
						sum_w = omt* omt* w[0] + 2 * omt * t * w[1] + t * t * w[2];
						
						x = (omt * omt * p0_r[0] * w[0] + 2 * omt * t * p1_r[0] * w[1] + t * t * p2_r[0] * w[2]) / sum_w;
						y = (omt * omt * p0_r[1] * w[0] + 2 * omt * t * p1_r[1] * w[1] + t * t * p2_r[1] * w[2]) / sum_w;

						t -= step;

						vertices.push_back(x);
						vertices.push_back(y);
					}
				}
				else {
					
					float L[3] = { A[1] - p1_r[1], p1_r[0] - A[0], (A[0] * p1_r[1] - p1_r[0] * A[1]) };
					
					float t1 = get_line_and_bezier_intersection_t(L, p0_r, p1_r, p2_r);
					

					float* p0_l = &u.second.second[next_ind][2];
					float* p1_l = &u.second.second[i][4];
					float* p2_l = &u.second.second[i][6];

					L[0] = A[1] - p1_l[1] ;
					L[1] = p1_l[0] - A[0];
					L[2] = (A[0] * p1_l[1] - p1_l[0] * A[1]);

					float t2 = get_line_and_bezier_intersection_t(L, p0_l, p1_l, p2_l);
					{
						const int STEPS_COUNT = 10;
						float step1 = (1.f - t1) / STEPS_COUNT;
						float step2 = (1.f - t2) / STEPS_COUNT;
						if (M > 2) {
							float t = t1;
							float x1 = (1 - t) * (1 - t) * p0_r[0] + 2 * (1 - t) * t * p1_r[0] + t * t * p2_r[0];
							float y1 = (1 - t) * (1 - t) * p0_r[1] + 2 * (1 - t) * t * p1_r[1] + t * t * p2_r[1];

							t = t2;
							float x2 = (1 - t) * (1 - t) * p0_l[0] + 2 * (1 - t) * t * p1_l[0] + t * t * p2_l[0];
							float y2 = (1 - t) * (1 - t) * p0_l[1] + 2 * (1 - t) * t * p1_l[1] + t * t * p2_l[1];

							vertices.push_back(A[0]);
							vertices.push_back(A[1]);
							vertices.push_back(x1);
							vertices.push_back(y1);
							vertices.push_back(x2);
							vertices.push_back(y2);
						}
						for (int i = 0; i < STEPS_COUNT; ++i) {
							float t = t1 + step1 * i;

							float x1 = (1 - t) * (1 - t) * p0_r[0] + 2 * (1 - t) * t * p1_r[0] + t * t * p2_r[0];
							float y1 = (1 - t) * (1 - t) * p0_r[1] + 2 * (1 - t) * t * p1_r[1] + t * t * p2_r[1];

							t = t2 + step2 * i;
							float x2 = (1 - t) * (1 - t) * p0_l[0] + 2 * (1 - t) * t * p1_l[0] + t * t * p2_l[0];
							float y2 = (1 - t) * (1 - t) * p0_l[1] + 2 * (1 - t) * t * p1_l[1] + t * t * p2_l[1];

							t = t1 + step1 * (i + 1);
							float x3 = (1 - t) * (1 - t) * p0_r[0] + 2 * (1 - t) * t * p1_r[0] + t * t * p2_r[0];
							float y3 = (1 - t) * (1 - t) * p0_r[1] + 2 * (1 - t) * t * p1_r[1] + t * t * p2_r[1];

							t = t2 + step2 * (i + 1);
							float x4 = (1 - t) * (1 - t) * p0_l[0] + 2 * (1 - t) * t * p1_l[0] + t * t * p2_l[0];
							float y4 = (1 - t) * (1 - t) * p0_l[1] + 2 * (1 - t) * t * p1_l[1] + t * t * p2_l[1];

							vertices.push_back(x1);
							vertices.push_back(y1);
							vertices.push_back(x3);
							vertices.push_back(y3);
							vertices.push_back(x4);
							vertices.push_back(y4);

							vertices.push_back(x1);
							vertices.push_back(y1);
							vertices.push_back(x4);
							vertices.push_back(y4);
							vertices.push_back(x2);
							vertices.push_back(y2);
						}
					}

				}
			}
		}
		
	}
#endif
}