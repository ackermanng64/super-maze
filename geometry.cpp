#include <cstdio>
#include <cmath>
#include "geometry.h"
void get_min_max_pos(float mmx[2], float mmy[2], float* pos, int N, float delta) {
	mmx[0] = mmx[1] = pos[0];
	mmy[0] = mmy[1] = pos[1];
	for (int i = 0; i < N; ++i) {
		if (mmx[0] > pos[2 * i]) {
			mmx[0] = pos[2 * i];
		}
		if (mmx[1] < pos[2 * i]) {
			mmx[1] = pos[2 * i];
		}
		if (mmy[0] > pos[2 * i + 1]) {
			mmy[0] = pos[2 * i + 1];
		}
		if (mmy[1] < pos[2 * i + 1]) {
			mmy[1] = pos[2 * i + 1];
		}
	}

	mmx[0] -= delta;
	mmx[1] += delta;
	mmy[0] -= delta;
	mmy[1] += delta;
}
float sqr_dist(float* a, float* b) {
	return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]);
}
float sqr_len(float* a) {
	return (a[0] * a[0] + a[1] * a[1]);
}
float dot_prd(float* a, float* b) {
	return a[0] * b[0] + a[1] * b[1];
}
float get_perp_vec(float* v, float* out) {
	out[0] = -v[1]; out[0] = v[0];
}
float get_perp_vec(float* p1, float* p2, float* out) {
	out[0] = p1[1] - p2[1]; out[1] = p2[0] - p1[0];
}
bool line_and_line_intersection(float* p1, float* p2, float* q1, float* q2, float* out) {
	float a = p1[0] * p2[1] - p1[1] * p2[0];
	float b = q1[0] * q2[1] - q1[1] * q2[0];
	float c = (p1[0] - p2[0]) * (q1[1] - q2[1]) - (p1[1] - p2[1]) * (q1[0] - q2[0]);
	if (c == 0) {
		return false;
	}
	out[0] = (a * (q1[0] - q2[0]) - b * (p1[0] - p2[0])) / c;
	out[1] = (a * (q1[1] - q2[1]) - b * (p1[1] - p2[1])) / c;
	return true;
}
bool line_seg_and_seg_intesection(float* p1, float* p2, float* q1, float* q2, float* out) {
	if (line_and_line_intersection(p1, p2, q1, q2, out)) {
		float s1, s2;
		if (p1[0] == p2[0]) {
			s1 = (out[1] - p1[1]) / (p2[1] - p1[1]);
		}
		else {
			s1 = (out[0] - p1[0]) / (p2[0] - p1[0]);
		}
		if (q1[0] == q2[0]) {
			s2 = (out[1] - q1[1]) / (q2[1] - q1[1]);
		}
		else {
			s2 = (out[0] - q1[0]) / (q2[0] - q1[0]);
		}
		if (s1 > 0 && s1 < 1 && s2 > 0 && s2 < 1) {
			return true;
		}
	}
	return false;
}

float dist_from_point_to_line(float* p1, float* p2, float* q) {
	float a = p2[0] - p1[0];
	float b = p2[1] - p1[1];
	return fabs(a * (p1[1] - q[1]) - b * (p1[0] - q[0])) / sqrt(a*a+b*b);
}

// negative if p lies to the left of ab
float point_orientation(float* a, float* b, float* p) {
	return (p[0] - a[0]) * (b[1] - a[1]) - (b[0] - a[0]) * (p[1] - a[1]);
}

// determine where the given point(p) lies in triangle
// return -1 if outside
// return 0/1/2 if on one of the sides, 0 - ab, 1-bc, 2-ca
// return 3 if inside
// th - threshold for how much error is allowed when deciding whether lies on a side
int point_location_test_v_triangle(float* a, float* b, float* c, float* p, float th) {
	float s1 = point_orientation(a, b, p);
	if (!(s1 > th || -s1 > th)) {
		float v1[2] = { b[0] - a[0], b[1] - a[1] };
		float v2[2] = { p[0] - a[0], p[1] - a[1] };
		float d = dot_prd(v1, v2);
		if (d > 0 && sqr_len(v1) > sqr_len(v2)) {
			return 0;
		}
		return -1;
	}

	float s2 = point_orientation(b, c, p);
	if (!(s2 > th || -s2 > th)) {
		float v1[2] = { c[0] - b[0], c[1] - b[1] };
		float v2[2] = { p[0] - b[0], p[1] - b[1] };
		float d = dot_prd(v1, v2);
		if (d > 0 && sqr_len(v1) > sqr_len(v2)) {
			return 1;
		}
		return -1;
	}

	if ((s1 > 0 && s2 < 0) || (s1 < 0 && s2 > 0)) {
		return -1;
	}
	
	s2 = point_orientation(c, a, p);
	if (!(s2 > th || -s2 > th)) {
		float v1[2] = { a[0] - c[0], a[1] - c[1] };
		float v2[2] = { p[0] - c[0], p[1] - c[1] };
		float d = dot_prd(v1, v2);
		if (d > 0 && sqr_len(v1) > sqr_len(v2)) {
			return 1;
		}
		return -1;
	}

	if ((s1 > 0 && s2 < 0) || (s1 < 0 && s2 > 0)) {
		return -1;
	}
	
	return 3;
}

void get_tri_circumcirlce(float* a, float* b, float* c, float* out) {
	float m1[2] = { (a[0] + b[0]) / 2, (a[1] + b[1]) / 2 };
	float m2[2] = { (b[0] + c[0]) / 2, (b[1] + c[1]) / 2 };
	float v1[2] = { m1[0] - (b[1] - a[1]), m1[1] + (b[0] - a[0]) };
	float v2[2] = { m2[0] - (c[1] - b[1]), m2[1] + (c[0] - b[0]) };
	line_and_line_intersection(m1, v1, m2, v2, out);
}
void get_tri_incenter(float* a, float* b, float* c, float* out) {
	float l1 = sqrt(sqr_dist(a, b));
	float l2 = sqrt(sqr_dist(b, c));
	float l3 = sqrt(sqr_dist(c, a));
	float ls = l1 + l2 + l3;
	out[0] = (l1 * c[0] + l2 * a[0] + l3 * b[0]) / ls;
	out[1] = (l1 * c[1] + l2 * a[1] + l3 * b[1]) / ls;
};	

void count_ray_cast_intersection(float* p1, float* p2, float* q, int* counter) {
	float max_x = p1[0];
	if (p1[0] < p2[0]) {
		max_x = p2[0];
	}
	float min_y = p1[1];
	float max_y = p2[1];
	if (p1[1] > p2[1]) {
		min_y = p2[1];
		max_y = p1[1];
	}
	if (q[1] > min_y && q[1] <= max_y && q[0] < max_x) {
		if (p1[1] != p2[1]) {
			if (p1[0] == p2[0]) {
				++(*counter);
			}
			else {
				float x = (q[1] - p1[1]) * (p2[0] - p1[0]) / (p2[1] - p1[1]) + p1[0];
				if (q[0] <= x) {
					++(*counter);
				}
			}
		}
	}
}
bool is_point_inside_polygon(float* points, int N, float* q) {
	int counter = 0;
	for (int i = 0; i < N - 1; ++i) {
		count_ray_cast_intersection(&points[2 * i], &points[2 * i + 2], q, &counter);
	}
	count_ray_cast_intersection(&points[2 * N - 2], &points[0], q, &counter);
	return (counter % 2 == 1);
}

int jarvis_march(float* points, int N) {
	if (N < 4) {
		return N;
	}
	int ind = 0;
	int M = 0;
	float point_on_hull[2] = { points[0], points[1] };
	for (int i = 1; i < N; ++i) {
		if (points[2 * i] < point_on_hull[0]) {
			point_on_hull[0] = points[2 * i];
			point_on_hull[1] = points[2 * i + 1];
			ind = i;
		}
	}

	while(true) {
		points[2 * ind] = points[2 * M];
		points[2 * ind + 1] = points[2 * M + 1];
		points[2 * M] = point_on_hull[0];
		points[2 * M + 1] = point_on_hull[1];
		++M;
		ind = M;
		for (int i = M; i < N; ++i) {
			if (point_orientation(point_on_hull, &points[2 * ind], &points[2 * i]) < 0) {
				ind = i;
			}
		}
		if(M == N || point_orientation(point_on_hull, &points[2*ind], &points[0])  < 0) {
			break;
		}
		point_on_hull[0] = points[2 * ind];
		point_on_hull[1] = points[2 * ind + 1];
	}
	return M;
}