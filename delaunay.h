#pragma once
#include <cmath>
#include <vector>

struct tri {
	int a[3];
};
struct sorted_tri {
	int a[3];
	sorted_tri() {
		set(-10, -20, -30);
	}
	sorted_tri(int a0, int a1, int a2) {
		set(a0, a1, a2);
	}
	void set(int a0, int a1, int a2) {

		a[0] = a0;
		a[1] = a1;
		a[2] = a2;
		if (a[0] > a[1]) {
			int t = a[0];
			a[0] = a[1];
			a[1] = t;
		}
		if (a[1] > a[2]) {
			int t = a[1];
			a[1] = a[2];
			a[2] = t;
		}
		if (a[0] > a[1]) {
			int t = a[0];
			a[0] = a[1];
			a[1] = t;
		}
	}
};
struct sorted_tri_comparator {
	bool operator() (const sorted_tri& t1, const sorted_tri& t2) const{
		if (t1.a[0] == t2.a[0]) {
			if (t1.a[1] == t2.a[1]) {
				return (t1.a[2] < t2.a[2]);
			}
			return (t1.a[1] < t2.a[1]);
		}
		return (t1.a[0] < t2.a[0]);
	}
};

struct sorted_tri_pair {
	static sorted_tri_comparator stc;
	sorted_tri st1;
	sorted_tri st2;
	sorted_tri_pair() {}
	sorted_tri_pair(const sorted_tri& _st1, const sorted_tri& _st2) {
		set(_st1, _st2);
	}
	void set(const sorted_tri& _st1, const sorted_tri& _st2) {
		st1 = _st1;
		st2 = _st2;
		if (stc(_st2, _st1)) {
			st1 = _st2;
			st2 = _st1;
		}
	}
};

struct sorted_tri_pair_comparator {
	bool operator() (const sorted_tri_pair& p1, const sorted_tri_pair& p2) const {
		if (!sorted_tri_pair::stc(p1.st1, p2.st1) && !sorted_tri_pair::stc(p2.st1, p1.st1)) {
			return sorted_tri_pair::stc(p1.st2, p2.st2);
		}
		return sorted_tri_pair::stc(p1.st1, p2.st1);
	}
};
struct delaunay_tri_edge {
	int a, b;
	delaunay_tri_edge() {
		a = b = -1;
	}
	delaunay_tri_edge(int _a, int _b) {
		set(_a, _b);
	}
	void set(int _a, int _b) {
		a = _a;
		b = _b;
		if (_a > _b) {
			a = _b;
			b = _a;
		}
	}
};
struct delaunay_edge_comparator {
	bool operator() (const delaunay_tri_edge& e1, const delaunay_tri_edge& e2) const {
		if (e1.a == e2.a) {
			return e1.b < e2.b;
		}
		return e1.a < e2.a;
	}
};
 
std::vector<tri> bowyer_watson_triangualtion(float* points, int N);
void t_voronoi(float* points, int N, std::vector<float>& vertices, float* border, int K);