#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>
#include <stack>
#include <set>
#include <unordered_set>
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

void t_voronoi(float* points, int N, std::vector<float>& vertices, float* border, int K) {
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
				printf("not inside: [%d] (%f, %f)\n", i, points[2 * i], points[2 * i + 1]);
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

	/*printf("ems: %ld\n", edge_map.size());
	for (auto& u : edge_map) {
		if (u.first.a > -1 && u.first.b > -1) 
		//	(u.second.a[1]  > -1 || u.second.a[1] == -4) && (u.second.a[2] > -1 || u.second.a[2] == -4)) 
		{
			printf("e: (%d %d) v1: (%d) v2: (%d)\n", u.first.a, u.first.b, u.second.a[1], u.second.a[2]);
		}
	}*/

	float min_len_threshold = 0;

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


	PointSelectionOption point_selection_option = PointSelectionOption::Circumcenter;
	std::map<sorted_tri, std::vector<float>, sorted_tri_comparator> valid_vor_verts;
	std::map<delaunay_tri_edge, std::pair<sorted_tri, sorted_tri>, delaunay_edge_comparator> valid_edges;
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
			if (sqr_dist(c1, c2) > min_len_threshold * min_len_threshold) 
			{
				//sorted_tri st1(u.first.a, u.first.b, c);
				//sorted_tri st2(u.first.a, u.first.b, d);
				//valid_edges.insert({ u.first, {st1, st2} });
				//valid_vor_verts.insert({ st1, {c1[0], c1[1]} });
				//valid_vor_verts.insert({ st2, {c2[0], c2[1]} });

				//printf("cd: %d %d\n", c, d);
				if (c > -1 && d > -1) 
				{
					bool b1 = is_point_inside_polygon(border, K, c1);
					bool b2 = is_point_inside_polygon(border, K, c2);
					if (b1 && b2) 
					{
						sorted_tri st1(u.first.a, u.first.b, c);
						sorted_tri st2(u.first.a, u.first.b, d);
						valid_edges.insert({ u.first, {st1, st2} });
						valid_vor_verts.insert({ st1, {c1[0], c1[1]} });
						valid_vor_verts.insert({ st2, {c2[0], c2[1]} });
					}
					else {
						has_invalid_vert_set.insert(u.first.a);
						has_invalid_vert_set.insert(u.first.b);
					}
				}
				else {
#if 1
					has_invalid_vert_set.insert(u.first.a);
					has_invalid_vert_set.insert(u.first.b);				
#endif
				}
			}
		}
	}
	
	printf("has invalid vert:\n");
	for (auto& u : has_invalid_vert_set) {
		printf("%d\n", u);
	}

	float min_len = 1000;
	int start_ind;
	std::vector<std::vector<int>> nbrdata(N);
	for (auto& u : valid_edges) {
		float* c1 = valid_vor_verts[u.second.first].data();
		float* c2 = valid_vor_verts[u.second.second].data();
		if (sqr_dist(c1, c2) < min_len) {
			min_len = sqr_dist(c1, c2);
		}
		
		if (has_invalid_vert_set.find(u.first.a) == has_invalid_vert_set.end() && has_invalid_vert_set.find(u.first.b) == has_invalid_vert_set.end()) {
			start_ind = u.first.a;
			nbrdata[u.first.a].push_back(u.first.b);
			nbrdata[u.first.b].push_back(u.first.a);
		}
		
	}
	printf("min len: %.9f\n", min_len);
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
				valid_edges.erase(dte);
				break;
			}
		}
	}


	float hfw = 0.0015;

	if(0)
	for (auto& u : valid_edges) {
		float* c1 = valid_vor_verts[u.second.first].data();
		float* c2 = valid_vor_verts[u.second.second].data();
		
		float v[2];
		float vd = std::sqrt(sqr_dist(c1, c2));
		if (vd == 0) continue;
		v[0] = hfw * (c2[1] - c1[1]) / vd;
		v[1] = hfw * (c1[0] - c2[0]) / vd;
		float v1[2] = { hfw / 2 * (c1[0] - c2[0]) / vd, hfw / 2 * (c1[1] - c2[1]) / vd };

		vertices.push_back(c1[0] + v[0]);
		vertices.push_back(c1[1] + v[1]);
		vertices.push_back(c2[0] + v[0]);
		vertices.push_back(c2[1] + v[1]);
		vertices.push_back(c2[0] - v[0]);
		vertices.push_back(c2[1] - v[1]);

		vertices.push_back(c2[0] - v[0]);
		vertices.push_back(c2[1] - v[1]);
		vertices.push_back(c1[0] - v[0]);
		vertices.push_back(c1[1] - v[1]);
		vertices.push_back(c1[0] + v[0]);
		vertices.push_back(c1[1] + v[1]);
	}
	
	std::map<sorted_tri, std::vector<sorted_tri>, sorted_tri_comparator> pmap;
	for (auto& u : valid_edges) {
		auto p1 = u.second.first;
		auto p2 = u.second.second;
		auto it = pmap.find(p1);
		if (it == pmap.end()) {
			pmap.insert(it, { p1, {u.second.second} });
		}
		else {
			it->second.push_back(u.second.second);
		}
		it = pmap.find(p2);
		if (it == pmap.end()) {
			pmap.insert(it, { p2, {u.second.first} });
		}
		else {
			it->second.push_back(u.second.first);
		}
	}
	if(1)
	for (auto& u : pmap) {
		float* center = valid_vor_verts[u.first].data();
		int N = u.second.size();
		std::vector<float> indices(N);
		std::vector<float> pos;
		for (int i = 0; i < N; ++i) {
			indices[i] = i;
			float* v = valid_vor_verts[u.second[i]].data();
			float dx = v[0] - center[0];
			float dy = v[1] - center[1];
			float d = sqrt(dx * dx + dy * dy);
			pos.push_back(dx / d + center[0]);
			pos.push_back(dy / d + center[1]);
		}
		
		// perform jarvis march 
		if (N > 3) {
			int ind = 0;
			int M = 0;
			float point_on_hull[2] = { pos[0], pos[1] };
			for (int i = 1; i < N; ++i) {
				if (pos[2 * i] < point_on_hull[0]) {
					point_on_hull[0] = points[2 * i];
					point_on_hull[1] = points[2 * i + 1];
					ind = i;
				}
			}

			while (true) {
				int j = indices[ind];
				indices[ind] = indices[M];
				indices[M] = j;
				pos[2 * ind] = pos[2 * M];
				pos[2 * ind + 1] = pos[2 * M + 1];
				pos[2 * M] = point_on_hull[0];
				pos[2 * M + 1] = point_on_hull[1];
				++M;
				ind = M;
				for (int i = M; i < N; ++i) {
					if (point_orientation(point_on_hull, &pos[2 * ind], &pos[2 * i]) < 0) {
						ind = i;
					}
				}
				if (M == N || point_orientation(point_on_hull, &pos[2 * ind], &pos[0]) < 0) {
					break;
				}
				point_on_hull[0] = pos[2 * ind];
				point_on_hull[1] = pos[2 * ind + 1];
			}

			if (M != N) {
				printf("jarvis march: smth is wrong ret: %d size: %d\n", M, N);
			}
		}

		std::vector<sorted_tri> tmp(N);
		for (int i = 0; i < N; ++i) {
			tmp[i] = u.second[indices[i]];
		}
		u.second = tmp;
	}

	
	auto get_bisector_vec_with_mag = [](float* p1, float* p2, float* p3, float hfw, float* out) {
		float incenter[2];
		get_tri_incenter(p1, p2, p3, incenter);
		float a = sqrt(sqr_dist(p1, p2));
		float b = sqrt(sqr_dist(p1, p3));
		float c = sqrt(sqr_dist(p2, p3));
		float s = (a + b + c) / 2;
		float in_r = sqrt((s - a) * (s - b) * (s - c) / s);
		float d = sqrt(sqr_dist(incenter, p1));
		float mag = hfw / in_r;
		out[0] = mag * (incenter[0] - p1[0]);
		out[1] = mag * (incenter[1] - p1[1]);
	};
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
		
	std::set<sorted_tri_pair, sorted_tri_pair_comparator> done_set;
	sorted_tri_pair stp;
	if(1)
	for (auto& u : pmap) {
		float* oo = valid_vor_verts[u.first].data();
		float* pp = oo;
		auto& nlist = u.second;
		int N = nlist.size();
		for (int i = 0; i < N; ++i) {
#if 0
			bool skip_edge = false;
			for (int j = 0; j < 2 && !skip_edge; ++j) {
				for (int k = j + 1; k < 3 && !skip_edge; ++k) {
					for (int l = 0; l < 2 && !skip_edge; ++l) {
						for (int m = l + 1; m < 3; ++m) {
							if (u.first.a[j] == nlist[i].a[l] && u.first.a[k] == nlist[i].a[m]) {
								if (has_invalid_vert_set.find(u.first.a[3 - j - k]) == has_invalid_vert_set.end() &&
									has_invalid_vert_set.find(nlist[i].a[3 - l - m]) != has_invalid_vert_set.end()) {
									printf("(%d, %d)\n", u.first.a[3 - j - k], nlist[i].a[3 - l - m]);
									skip_edge = true;
									break;
								}
							}
						}
					}
				}
			}
			if (skip_edge) continue;
#endif

			stp.set(u.first, nlist[i]);
			auto it = done_set.find(stp);
			if (it != done_set.end()) {
				continue;
			}
			done_set.insert(it, stp);

			float* tt = valid_vor_verts[nlist[i]].data();
			float dx = pp[0] - tt[0];
			float dy = pp[1] - tt[1];
			if(dx*dx + dy*dy > 0) 
			{
				float u1[2] = {0, 0};
				float u2[2] = {0, 0};

				if (N == 1) {
					get_vs(pp, tt, nullptr, nullptr, hfw, u1, u2);
				}
				else if (N == 2) {
					get_vs(pp, tt, valid_vor_verts[nlist[(i - 1 + N) % N]].data(), nullptr, hfw, u1, u2);
				}
				else{
					get_vs(pp, tt, valid_vor_verts[nlist[(i - 1 + N) % N]].data(), valid_vor_verts[nlist[(i + 1) % N]].data(), hfw, u1, u2);
				}

				auto& mlist = pmap[nlist[i]];
				int M = mlist.size();
				for (int j = 0; j < M; ++j) {
					if (u.first.a[0] == mlist[j].a[0] && u.first.a[1] == mlist[j].a[1] && u.first.a[2] == mlist[j].a[2]) {
						float v1[2] = {0, 0};
						float v2[2] = {0, 0};

						if (M == 1) {
							get_vs(tt, pp, nullptr, nullptr, hfw, v1, v2);
						}
						else if (M == 2) {
							get_vs(tt, pp, valid_vor_verts[mlist[(j - 1 + M) % M]].data(), nullptr, hfw, v1, v2);
						}
						else {
							get_vs(tt, pp, valid_vor_verts[mlist[(j - 1 + M) % M]].data(), valid_vor_verts[mlist[(j + 1) % M]].data(), hfw, v1, v2);
						}

						float l1[2] = { pp[0] + u1[0], pp[1] + u1[1] };
						float l2[2] = { tt[0] + v1[0], tt[1] + v1[1] };
						float* w1 = v1, * w2 = v2;
						float res[2];
						if (line_seg_and_seg_intesection(l1, l2, pp, tt, res)) {
							w1 = v2;
							w2 = v1;
						}
						
						//printf("(%.3f, %.3f), (%.3f, %.3f), (%.3f, %.3f), (%.3f, %.3f)\n", 
						//	1e4 * u1[0], 1e4 * u1[1], 1e4 * u2[0], 1e4 * u2[1], 1e4 * w1[0], 1e4 * w1[1], 1e4 * w2[0], 1e4* w2[1]);
						if (1) {
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
						}else {
							/*vertices.push_back(pp.first);
							vertices.push_back(pp.second);
							vertices.push_back(tt.first);
							vertices.push_back(tt.second);
							
							vertices.push_back(pp.first);
							vertices.push_back(pp.second);
							vertices.push_back(pp.first + u1[0]);
							vertices.push_back(pp.second + u1[1]);

							vertices.push_back(tt.first);
							vertices.push_back(tt.second);
							vertices.push_back(tt.first + w1[0]);
							vertices.push_back(tt.second + w1[1]);
							
							vertices.push_back(pp.first);
							vertices.push_back(pp.second);
							vertices.push_back(pp.first + u2[0]);
							vertices.push_back(pp.second + u2[1]);

							vertices.push_back(tt.first);
							vertices.push_back(tt.second);
							vertices.push_back(tt.first + w2[0]);
							vertices.push_back(tt.second + w2[1]);*/
						}
						break;
					}
				}
			}
		}
	}
}