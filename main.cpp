#include <SFML/Graphics.hpp>
#include <vector>
#include <chrono>
#include <queue>
#include "geometry.h"
#include "delaunay.h"
#include "voronoi.h"
#include "low_disc_seq.h"
using namespace std;

void rescale_pos(float mmx[2], float mmy[2], float* pos, int N, float width, float height, float hmargin, float vmargin) {
    float max_size = fmaxf(mmx[1] - mmx[0], mmy[1] - mmy[0]);
    for (int i = 0; i < N; ++i) {
        pos[2 * i] = (width - 2*hmargin) * (pos[2*i] - mmx[0]) / (mmx[1] - mmx[0]) + hmargin;
        pos[2 * i + 1] = (height - 2*vmargin) * (pos[2*i + 1] - mmy[0]) / (mmy[1] - mmy[0]) + vmargin;
        pos[2 * i + 1] = height - pos[2 * i + 1];
    }
}
void rescale_pos(float mmx[2], float mmy[2], vector<sf::Vertex>& vs, float width, float height, float hmargin, float vmargin) {
    float max_size = fmaxf(mmx[1] - mmx[0], mmy[1] - mmy[0]);
    for (int i = 0; i < vs.size(); ++i) {
        vs[i].position.x = (width - 2*hmargin) * (vs[i].position.x - mmx[0]) / (mmx[1] - mmx[0]) + hmargin;
        vs[i].position.y = (height - 2*vmargin) * (vs[i].position.y - mmy[0]) / (mmy[1] - mmy[0]) + vmargin;
        vs[i].position.y = height - vs[i].position.y;
    }
}

void test_point_location_test_v_triangle() {
    float t[] = { -1.9062, -1.6528, 2.5352, -1.3018, 2.2206, 2.2683 };
    float p1[] = { -0.7, -1.17 };
    float p2[] = { 1.81, -1.36 };
    float p3[] = { 1.38, -1.33 };

    printf("p1 %d\n", point_location_test_v_triangle(&t[0], &t[2], &t[4], p1, 0.01));
    printf("p2 %d\n", point_location_test_v_triangle(&t[0], &t[2], &t[4], p2, 0.01));
    printf("p3 %d\n", point_location_test_v_triangle(&t[0], &t[2], &t[4], p3, 0.01));
}

void extract_triangles(float* points, vector<tri>& triangles, vector<sf::Vertex>& vertices, bool print) {
    
    for (auto& u : triangles) {
        if (print) {
            printf("tri: (%d, %d, %d)\n", u.a[0], u.a[1], u.a[2]);
        }
        
        sf::Vertex v;
        v.color = sf::Color(rand() % 256, rand() % 256, rand() % 256);
        //v.color = sf::Color::White;
        v.position = sf::Vector2f(points[2*u.a[0]], points[2 * u.a[0] + 1]);
        vertices.push_back(v);
        v.position = sf::Vector2f(points[2 * u.a[1]], points[2 * u.a[1] + 1]);
        vertices.push_back(v);

        v.position = sf::Vector2f(points[2 * u.a[1]], points[2 * u.a[1] + 1]);
        vertices.push_back(v);
        v.position = sf::Vector2f(points[2 * u.a[2]], points[2 * u.a[2] + 1]);
        vertices.push_back(v);

        v.position = sf::Vector2f(points[2 * u.a[2]], points[2 * u.a[2] + 1]);
        vertices.push_back(v);
        v.position = sf::Vector2f(points[2 * u.a[0]], points[2 * u.a[0] + 1]);
        vertices.push_back(v);
        
    }
}
vector<float> test_delaunay1(vector<sf::Vertex>& triangle_vertices) {
    vector<float> sites = {
       0, 0,
       0.18, 0.76,
       0.75, 0.18,
       0.62, 0.62
    };
    
    auto T = bowyer_watson_triangualtion(sites.data(), sites.size() / 2);
    extract_triangles(sites.data(), T, triangle_vertices, true);
    return sites;
}

vector<float> test_delaunay2(vector<sf::Vertex>& triangle_vertices) {
    vector<float> sites = {
       4.5, 5.21,
       3, 3.5,
       3.89, 3.73,
       4.94, 2.82,
       4.37, 4.48,
       3.73, 3.08,
       3.5, 3.5,
       4.28, 3.45,
       3.09, 2.99,
       4.13, 2.76,
       3.46, 3.24,
       3.88, 3.41
    };

    auto T = bowyer_watson_triangualtion(sites.data(), sites.size() / 2);
    extract_triangles(sites.data(), T, triangle_vertices, true);
    return sites;
}

vector<float> test_delaunay3(vector<sf::Vertex>& triangle_vertices) {
    vector<float> sites;

    if(0)
    for (int i = 0; i < 200; ++i) {
        sites.push_back((rand() % 2560) / 2559.f);
        sites.push_back((rand() % 2560) / 2559.f);
    }

    auto pg = poission_disc_sampling(1, 1, 0.02);
    if(1)
    for (auto& u : pg) {
        for (auto& v : u) {
            if (v.first != -1) {
                sites.push_back(v.first);
                sites.push_back(v.second);
            }
        }
    }
    if(0)
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            float x = i / 10.f + 0.03f;
            float y = j / 10.f + 0.03;
            sites.push_back(x);
            sites.push_back(y);
        }
    }
    if(0) sites = {
        0.28, 0.61,
        0.33, 0.39,
        0.58, 0.38,
    };
    auto start = std::chrono::steady_clock::now();
    auto T = bowyer_watson_triangualtion(sites.data(), sites.size() / 2);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    printf("bw1: %lfms\n", elapsed.count());
    //auto T = t_voronoi(sites.data(), sites.size() / 2);
    extract_triangles(sites.data(), T, triangle_vertices, false);
    return sites;
}

void test_voronoi1(vector<float>& sites, vector<sf::Vertex>& edges, vector<float>& border,
    std::map<sorted_tri, std::vector<float>, sorted_tri_comparator>& valid_vor_verts,
    std::vector<std::vector<sorted_tri_pair>>& site_polygon,
    std::vector<std::vector<int>>& nbrdata,
    std::map<delaunay_tri_edge, std::pair<sorted_tri, sorted_tri>, delaunay_edge_comparator>& valid_edge,
    std::set<sorted_tri_pair, sorted_tri_pair_comparator>& removed_n_edges,
    int& start_ind,
    float start_pos[2]) { 
    vector<float> vertices;
    auto start = std::chrono::steady_clock::now();
    t_voronoi(  sites.data(), sites.size() / 2, vertices, border.data(), border.size() / 2, valid_vor_verts, site_polygon, nbrdata,
                valid_edge, removed_n_edges, start_ind, start_pos);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    printf("bw2: %lfms\n", elapsed.count());
    printf("site count: %lld\n", (sites.size() / 2));
    sf::Vertex v;
    v.color = sf::Color::White;
    for (int i = 0; i < vertices.size(); i += 2) {
        v.position.x = vertices[i];
        v.position.y = vertices[i + 1];
        edges.push_back(v);
    }
    printf("triangles count: %lld\n", (edges.size() / 3));
}
struct st_pn {
    sorted_tri prev;
    sorted_tri next;
    st_pn() {
        prev = next = sorted_tri();
    }
    st_pn(sorted_tri p, sorted_tri n) {
        prev = p;
        next = n;
    }
};
void f(int* p) {
    int* q = p;
    printf("%p %p\n", p, q);
}
void t_voronoi_post_process(std::map<sorted_tri, std::vector<float>, sorted_tri_comparator>& valid_vor_verts, float start_pos[2],
                            std::vector<std::vector<sorted_tri_pair>>& site_polygon, 
                            std::vector<std::vector<sorted_tri>>& site_sorted_polygon,
                            float mmx[2], float mmy[2], float width, float height, float hmargin, float vmargin) {

    float max_size = fmaxf(mmx[1] - mmx[0], mmy[1] - mmy[0]);
    for (auto& u : valid_vor_verts) {
        u.second[0] = (width - 2 * hmargin) * (u.second[0] - mmx[0]) / (mmx[1] - mmx[0]) + hmargin;
        u.second[1] = (height - 2 * vmargin) * (u.second[1] - mmy[0]) / (mmy[1] - mmy[0]) + vmargin;
        u.second[1] = height - u.second[1];
    }

    start_pos[0] = (width - 2 * hmargin) * (start_pos[0]  - mmx[0]) / (mmx[1] - mmx[0]) + hmargin;
    start_pos[1] = (height - 2 * vmargin) * (start_pos[1] - mmy[0]) / (mmy[1] - mmy[0]) + vmargin;
    start_pos[1] = height - start_pos[1];

    // sort edges
   
    site_sorted_polygon.resize(site_polygon.size());
    for (int i = 0; i < site_polygon.size(); ++i) {        
        auto u = site_polygon[i];

        if (u.size() == 0) continue;
        int size = u.size();

        auto& v = site_sorted_polygon[i];
        std::map<sorted_tri, std::vector<sorted_tri>, sorted_tri_comparator> nmap;
        std::set<sorted_tri, sorted_tri_comparator> visited_set;
        for (auto& v : u) {
            nmap[v.st1].push_back(v.st2);
            nmap[v.st2].push_back(v.st1);
        }

        auto next = u[0].st1;
        visited_set.insert(next);
        v.push_back(next);

        while (visited_set.size() != u.size()) {
            auto& arr = nmap[next];
            if (visited_set.find(arr[0]) == visited_set.end()) {
                visited_set.insert(arr[0]);
                next = arr[0];

            }
            else if(visited_set.find(arr[1]) == visited_set.end()) {
                visited_set.insert(arr[1]);
                next = arr[1];
            }
            v.push_back(next);
        }
    }
}
struct player_location_data {
    int seg_ind = 0;
    float t;
    std::vector<float> crossing_ts;
    std::vector<int> others;
};

int main()
{
    std::vector<int> seeds_with_bugs = {
        1606805185,
        1606805427,
        1607858997,
        1607859289,
        1607944581, // border cull bug with circumcenter
        1608274936, // render problem -almost solved
        1608464112, // render problem2
    };
    time_t seed = time(NULL);
    srand(1608464112);
    printf("seed: %ld\n", seed);
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    int window_width = 1000;
    int window_height = 1000;
    int hmargin = 10;
    int vmargin = 10;
    sf::RenderWindow window(sf::VideoMode(window_width, window_height), "voronoi", sf::Style::Default, settings);
    sf::CircleShape shape(5);
    vector<sf::Vertex> triangles_vertices;
    vector<sf::Vertex> voronoi_edges;
    
    vector<float> cull_polygon1 = { 0.3, 0.3, 0.6, 0.6, 0.45, 0.75 };
    //vector<float> cull_polygn2 = { 0.3, 0.3, 0.6, 0.6, 0.45, 0.75, 0.12, 0.75, 0.07, 0.56, 0.24, 0.49 };
    vector<float> cull_polygon5 = 
    { 
        0.2, 0.2, 
        0.28, 0.26, 
        0.32, 0.4, 
        0.5, 0.3, 
        0.47, 0.22, 
        0.58, 0.21,
        0.58, 0.64,
        0.18, 0.66
    };
    vector<float> cull_polygon3 =
    {
        0.02, 0.02,
        0.98, 0.02,
        0.98, 0.98,
        0.02, 0.98
    };
    vector<float> cull_polygon =
    {
        0, 0,
        1, 0,
        1, 1,
        0, 1
    };
    vector<float> cull_polygon2 =
    {
        0.12, 0.11,
        0.39, 0.07,
        0.48, 0.29,
        0.61, 0.52,
        0.59, 0.74,
        0.44, 0.89,
        0.22, 0.89,
        0.34, 0.75,
        0.44, 0.62,
        0.27, 0.57,
        0.12, 0.61,
        0.09, 0.45,
        0.23, 0.4,
        0.34, 0.31,
        0.24, 0.2
    };
    vector<sf::Vertex> cull_polygon_verts; 
    {
        int S = cull_polygon.size() / 2;
        for (int i = 0; i < S; ++i) {
            sf::Vertex v;
            v.color = sf::Color::Red;
            v.position.x = cull_polygon[2 * i];
            v.position.y = cull_polygon[2 * i + 1];
            cull_polygon_verts.push_back(v);
            v.position.x = cull_polygon[2 * ((i + 1) % S)];
            v.position.y = cull_polygon[2 * ((i + 1) % S) + 1];
            cull_polygon_verts.push_back(v);
        }
    }
    // sites are generted in [0,1]x[0,1] grid
    vector<float> sites;
    sites = test_delaunay3(triangles_vertices);
    vector<sf::Color> site_colors(sites.size() / 2);
    for (auto& u : site_colors) {
        u = sf::Color(rand() % 256, rand() % 256, rand() % 256);
    }

    std::map<sorted_tri, std::vector<float>, sorted_tri_comparator> valid_vor_verts;
    std::vector<std::vector<sorted_tri_pair>> site_polygon;
    std::vector<std::vector<int>> nbrdata;
    std::map<delaunay_tri_edge, std::pair<sorted_tri, sorted_tri>, delaunay_edge_comparator> valid_edges;
    std::set<sorted_tri_pair, sorted_tri_pair_comparator> removed_edges_set;
    std::vector<std::vector<sorted_tri>> site_sorted_polygon;
    int start_ind;
    float start_pos[2];
    test_voronoi1(sites, voronoi_edges, cull_polygon, valid_vor_verts, site_polygon, nbrdata, valid_edges, removed_edges_set, start_ind, start_pos);
    float mmx[2], mmy[2];
    get_min_max_pos(mmx, mmy, sites.data(), sites.size() / 2, 0.075f);
    t_voronoi_post_process(valid_vor_verts, start_pos, site_polygon, site_sorted_polygon, mmx, mmy, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, triangles_vertices, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, sites.data(), sites.size() / 2, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, voronoi_edges, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, cull_polygon.data(), cull_polygon.size() / 2, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, cull_polygon_verts, window_width, window_height, hmargin, vmargin);

    printf("start_pos: (%f, %f)\n", start_pos[0], start_pos[1]);
    printf("%ld\n", voronoi_edges.size() / 2);

    sf::Font font;
    if (!font.loadFromFile("Roboto-Bold.ttf")) {
        printf("could not load font file\n");
    }
    sf::Text text;
    text.setFont(font);
    text.setFillColor(sf::Color::Red);
    text.setCharacterSize(12);

    bool display_del = true;
    bool display_vor = true;
    bool display_sites = true;
    bool display_cull_polygon = true;

    sf::ConvexShape cs;

    std::vector<sf::Vertex> start_ind_polygon_data;
    
    for (auto& u : site_polygon[start_ind]) {
        sf::Vertex v;
        v.color = sf::Color::Red;
        float* p1 = valid_vor_verts[u.st1].data();
        v.position.x = p1[0];
        v.position.y = p1[1];
        start_ind_polygon_data.push_back(v);
     
        float* p2 = valid_vor_verts[u.st2].data();
        v.position.x = p2[0];
        v.position.y = p2[1];
        start_ind_polygon_data.push_back(v);

        //printf("(%d, %d, %d) (%d, %d, %d)\n", u.st1.a[0], u.st1.a[1], u.st1.a[2], u.st2.a[0], u.st2.a[1], u.st2.a[2]);
    }

    float area = 0;
    float centroid[2] = { 0, 0 };
    float cx = 0;
    float cy = 0;

    for (int i = 0; i < site_sorted_polygon[start_ind].size(); ++i) {
        //printf("(%d, %d, %d)\n", site_sorted_polygon[start_ind][i].a[0], site_sorted_polygon[start_ind][i].a[1], site_sorted_polygon[start_ind][i].a[2]);
        int j = (i + 1) % site_sorted_polygon[start_ind].size();
        float* p1 = valid_vor_verts[site_sorted_polygon[start_ind][i]].data();
        float* p2 = valid_vor_verts[site_sorted_polygon[start_ind][j]].data();

        float s = p1[0] * p2[1] - p1[1] * p2[0];
        area += s;
        cx += (p1[0] + p2[0]) * s;
        cy += (p1[1] + p2[1]) * s;
    }
    area /= 2;
    cx /= (6 * area);
    cy /= (6 * area);
    centroid[0] = cx;
    centroid[1] = cy;

    player_location_data pld;
    bool is_mouse_left_held = false;
    float mouse_press_pos[2];
    float displacement_vec[2];
    std::vector<sf::Vertex> movement_path_line;
    std::vector<sf::Vertex> path_to_nbrs_data;
    
    const int max_depth = 3;
    std::map<int, int> bfs_prev_map;
    std::queue<int> bfs_queue;

    for (auto& u : nbrdata[start_ind]) {
        printf("%d\n", u);
    }

    auto set_path_to_nbrs = [&]() {
        int depth = 0;
        bfs_prev_map.clear();
        while (!bfs_queue.empty()) {
            bfs_queue.pop();
        }

        pld.seg_ind = -1;
        pld.crossing_ts.clear();
        pld.others.clear();

        path_to_nbrs_data.clear();

        bfs_prev_map.insert({ start_ind, -1 });
        bfs_queue.push(start_ind);
        int curr_count = 1;
        int next_count = 0; 

        while (!bfs_queue.empty()) {
            int ind = bfs_queue.front();
            bfs_queue.pop();
            --curr_count;

            for (int i = 0; i < nbrdata[ind].size(); ++i) {
                int other = nbrdata[ind][i];
                auto& pst = valid_edges[delaunay_tri_edge(ind, other)];
                sorted_tri_pair stp(pst.first, pst.second);
                sorted_tri& st1 = stp.st1;
                sorted_tri& st2 = stp.st2;
                float* p1 = valid_vor_verts[st1].data();
                float* p2 = valid_vor_verts[st2].data();
                if (removed_edges_set.find(stp) != removed_edges_set.end()) {
                    if (bfs_prev_map.find(other) == bfs_prev_map.end()) {
                        bfs_queue.push(other);
                        bfs_prev_map.insert({ other, ind });
                        ++next_count;
                        
                        if (ind == start_ind) {
                            float* q1 = &sites[2 * ind];
                            float* q2 = &sites[2 * other];
                            float res[2];
                            if (line_seg_and_seg_intesection(p1, p2, q1, q2, res)) {
                                if (fabs(q1[0] - q2[0]) > fabs(q1[1] - q2[1])) {
                                    pld.crossing_ts.push_back((res[0] - q1[0]) / (q2[0] - q1[0]));
                                }
                                else {
                                    pld.crossing_ts.push_back((res[1] - q1[1]) / (q2[1] - q1[1]));
                                }
                            }
                            pld.others.push_back(other);
                        }

                        sf::Vertex v;
                        v.color = sf::Color::Red;

                        v.position.x = sites[2 * ind];
                        v.position.y = sites[2 * ind + 1];
                        path_to_nbrs_data.push_back(v);
                        v.position.x = sites[2 * other];
                        v.position.y = sites[2 * other + 1];
                        path_to_nbrs_data.push_back(v);

                    }
                }
            }
            
            if (curr_count == 0) {
                ++depth;
                curr_count = next_count;
                next_count = 0;
            }
            if (depth == max_depth) {
                break;
            }
        }
    };
    set_path_to_nbrs();

    int move_target;
    float move_t;
    bool is_target_reached = true;

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::A) {
                    display_del = !display_del;
                }
                if (event.key.code == sf::Keyboard::S) {
                    display_vor = !display_vor;
                }
                if (event.key.code == sf::Keyboard::D) {
                    display_sites = !display_sites;
                }
                if (event.key.code == sf::Keyboard::W) {
                    display_cull_polygon = !display_cull_polygon;
                }
            }
            if (event.type == sf::Event::MouseButtonPressed && is_target_reached) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    is_mouse_left_held = true;
                    mouse_press_pos[0] = event.mouseButton.x;
                    mouse_press_pos[1] = event.mouseButton.y;
                    float x = mouse_press_pos[0];
                    float y = mouse_press_pos[1];

                    for (auto& u : bfs_prev_map) {
                        if (u.first == start_ind) continue;

                        int counter = 0;
                        int i;
                        float xinters;
                        
                        for (int i = 0, j; i < site_sorted_polygon[u.first].size(); ++i) {
                            j = (i + 1 == site_sorted_polygon[u.first].size() ? 0 : i + 1);
                            sorted_tri& st1 = site_sorted_polygon[u.first][i];
                            sorted_tri& st2 = site_sorted_polygon[u.first][j];

                            float* p1 = valid_vor_verts[st1].data();
                            float* p2 = valid_vor_verts[st2].data();

                            if (y > fmin(p1[1], p2[1]) && y <= fmax(p1[1], p2[1]) && x <= fmax(p1[0], p2[0])) {
                                if (p1[1] != p2[1]) {
                                    xinters = (y - p1[1]) * (p2[0] - p1[0]) / (p2[1] - p1[1]) + p1[0];
                                    if (p1[0] == p2[0] || x <= xinters) {
                                        ++counter;
                                    }
                                }
                            }
                        }
                        if (counter % 2 == 1) {
                            //start_ind = u.first;
                            //start_pos[0] = sites[2 * start_ind];
                            //start_pos[1] = sites[2 * start_ind + 1];
                            //set_path_to_nbrs();
                            move_t = 0;
                            move_target = u.first;
                            is_target_reached = false;
                            
                            int prev = u.first;
                            int next = bfs_prev_map[prev];
                            while (next != start_ind) {
                                prev = next;
                                next = bfs_prev_map[next];
                            }
                            
                            for (int i = 0; i < pld.others.size(); ++i) {
                                if (prev == pld.others[i]) {
                                    pld.seg_ind = i;
                                    break;
                                }
                            }
                            
                            break;
                        }
                    }
                    
                }
            }
            if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    is_mouse_left_held = false;
                }
            }
            if (event.type == sf::Event::MouseMoved && is_mouse_left_held)
            {
                
            }
        }

        if (!is_target_reached) {
            if (move_t >= 1) {
                start_ind = pld.others[pld.seg_ind];
                set_path_to_nbrs();
                if (start_ind == move_target) {
                    is_target_reached = true;
                }
                else {
                    int prev = move_target;
                    int next = bfs_prev_map[prev];
                    while (next != start_ind) {
                        prev = next;
                        next = bfs_prev_map[next];
                    }

                    for (int i = 0; i < pld.others.size(); ++i) {
                        if (prev == pld.others[i]) {
                            pld.seg_ind = i;
                            break;
                        }
                    }
                    move_t = 0;
                }
            }
            else {
                move_t += 0.005;
                int other = pld.others[pld.seg_ind];
                start_pos[0] = (sites[2 * other] - sites[2 * start_ind]) * move_t + sites[2 * start_ind];
                start_pos[1] = (sites[2 * other + 1] - sites[2 * start_ind + 1]) * move_t + sites[2 * start_ind + 1];
            }
        }

        window.clear();
        if (display_del) {
            window.draw(triangles_vertices.data(), triangles_vertices.size(), sf::Lines);
        }
        if (display_vor) {
            window.draw(voronoi_edges.data(), voronoi_edges.size(), sf::Triangles);
        }
        if (display_cull_polygon) {
            window.draw(cull_polygon_verts.data(), cull_polygon_verts.size(), sf::Lines);
        }
        for (int i = 0; i < sites.size() && display_sites; i += 2) {
            float x = sites[i];
            float y = sites[i + 1];
            shape.setPosition(sf::Vector2f(x - 5, y - 5));
            shape.setFillColor(site_colors[i/2]);
            window.draw(shape);
            //if (bfs_prev_map.find(i / 2) != bfs_prev_map.end()) 
            {
                text.setPosition(x, y - 20);
                text.setString(to_string(i / 2));
                window.draw(text);
            }
        }
        for (int i = 0; i < cull_polygon.size() && display_cull_polygon; i += 2) {
            float x = cull_polygon[i];
            float y = cull_polygon[i + 1];
            shape.setPosition(sf::Vector2f(x - 5, y - 5));
            shape.setFillColor(sf::Color::Blue);
            window.draw(shape);
            text.setPosition(x, y - 20);
            text.setString("B" + to_string(i / 2));
            window.draw(text);
        }

        
        {
            shape.setPosition(sf::Vector2f(start_pos[0] - 5, start_pos[1] - 5));
            shape.setFillColor(sf::Color::Red);
            window.draw(shape);
        }
        {
            window.draw(movement_path_line.data(), movement_path_line.size(), sf::Lines);
        }
        {
            window.draw(path_to_nbrs_data.data(), path_to_nbrs_data.size(), sf::Lines);
        }
        
        window.display();
    }
    return 0;
}