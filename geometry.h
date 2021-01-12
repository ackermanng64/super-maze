#pragma once
void get_min_max_pos(float mmx[2], float mmy[2], float* pos, int N, float delta);
float sqr_dist(float* a, float* b);
inline float sqr_len(float* a);
inline float dot_prd(float* a, float* b);
// perp vector with the end to the left of v
inline float get_perp_vec(float* v, float* out);
// perp vector with the end to the left of p1p2
inline float get_perp_vec(float* p1, float* p2, float* out);
// line1: defined by p1, p2; line2: defined by q1, q2; make sure lines are not parallel - no check is done inside the function
bool line_and_line_intersection(float* p1, float* p2, float* q1, float* q2, float* out);
bool line_seg_and_seg_intesection(float* p1, float* p2, float* q1, float* q2, float* out);
bool line_and_line_intersection(int* p1, int* p2, int* q1, int* q2, float* out);
bool line_seg_and_seg_intesection(int* p1, int* p2, int* q1, int* q2, float* out);
// sqr distance from point q to a line defined with points p1, p2 
float sqr_dist_from_point_to_line(float* p1, float* p2, float* q);
// distance from point q to a line defined with points p1, p2 
float dist_from_point_to_line(float* p1, float* p2, float* q);
inline float point_orientation(float* a, float* b, float* p);
inline int point_orientation(int* a, int* b, int* p);
int point_location_test_v_triangle(float* a, float* b, float* c, float* p, float th);
void get_tri_circumcirlce(float* a, float* b, float* c, float* out);
void get_tri_incenter(float* a, float* b, float* c, float* out);
void count_ray_cast_intersection(float* p1, float* p2, float* q, int* counter);
bool is_point_inside_polygon(float* points, int N, float* q);
// points get mutated in this method, first M (return of the function) are convex hull
int jarvis_march(float* points, int N);
// sort vectors vector in CCW direction
void sort_vec_by_angle(float* vs, int N);