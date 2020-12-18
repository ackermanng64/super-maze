#include <SFML/Graphics.hpp>
#include <vector>
#include <chrono>
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

    auto pg = poission_disc_sampling(1, 1, 0.035);
    if(1)
    for (auto& u : pg) {
        for (auto& v : u) {
            if (v.first != -1) {
                sites.push_back(v.first);
                sites.push_back(v.second);
            }
        }
    }
    /*for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            float x = i / 10.f;
            float y = j / 10.f;
            sites.push_back(x);
            sites.push_back(y);
        }
    }*/
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

void test_voronoi1(vector<float>& sites, vector<sf::Vertex>& edges, vector<float>& border) {
    vector<float> vertices;
    auto start = std::chrono::steady_clock::now();
    t_voronoi(sites.data(), sites.size() / 2, vertices, border.data(), border.size() / 2);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    printf("bw2: %lfms\n", elapsed.count());

    sf::Vertex v;
    v.color = sf::Color::White;
    for (int i = 0; i < vertices.size(); i += 2) {
        v.position.x = vertices[i];
        v.position.y = vertices[i + 1];
        edges.push_back(v);
    }
    printf("triangles count: %lld\n", (edges.size() / 3));
}

int main()
{
    std::vector<int> seeds_with_bugs = {
        1606805185,
        1606805427,
        1607858997,
        1607859289,
        1607944581 // border cull bug with circumcenter
    };
    time_t seed = time(NULL);
    srand(1608274936);
    printf("seed: %ld\n", seed);
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    int window_width = 900;
    int window_height = 900;
    int hmargin = 100;
    int vmargin = 100;
    sf::RenderWindow window(sf::VideoMode(window_width, window_height), "voronoi", sf::Style::Default, settings);
    sf::CircleShape shape(5);
    vector<sf::Vertex> triangles_vertices;
    vector<sf::Vertex> voronoi_edges;
    
    vector<float> cull_polygon1 = { 0.3, 0.3, 0.6, 0.6, 0.45, 0.75 };
    //vector<float> cull_polygon2 = { 0.3, 0.3, 0.6, 0.6, 0.45, 0.75, 0.12, 0.75, 0.07, 0.56, 0.24, 0.49 };
    vector<float> cull_polygon3 = 
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
    vector<float> cull_polygon4 =
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

    test_voronoi1(sites, voronoi_edges, cull_polygon);
    float mmx[2], mmy[2];
    get_min_max_pos(mmx, mmy, sites.data(), sites.size() / 2, 0.075f);
    rescale_pos(mmx, mmy, triangles_vertices, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, sites.data(), sites.size() / 2, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, voronoi_edges, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, cull_polygon.data(), cull_polygon.size() / 2, window_width, window_height, hmargin, vmargin);
    rescale_pos(mmx, mmy, cull_polygon_verts, window_width, window_height, hmargin, vmargin);


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
            text.setPosition(x, y - 20);
            text.setString(to_string(i/2));
            window.draw(text);
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
        window.display();
    }

    return 0;
}