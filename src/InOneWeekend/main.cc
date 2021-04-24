//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

// needed for tiny obj loader
#define TINYOBJLOADER_IMPLEMENTATION

#include "rtweekend.h"

#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include "triangle.h"

#include <iostream>

#include <vector>
#include <string>
#include "external/tinyobjloader/tiny_obj_loader.h"

bool read_object_file(
    std::string obj_fname, 
    hittable_list& world, 
    const double scale = 1.0, 
    const shared_ptr<material> obj_material = nullptr, 
    const vec3 translation_vector = vec3(0,0,0)
    ) {
    // if (world == NULL) world = new hittable_list();

    tinyobj::ObjReaderConfig reader_config; //triangulate is automatically set to true
    tinyobj::ObjReader reader;
    // Spit out errors/warnings if you get some, terminate if error
    if (!reader.ParseFromFile(obj_fname, reader_config)) {
        if (!reader.Error().empty()) {std::cerr << "TinyObjReader Error: " << reader.Error() << std::endl;}
        return false;
    }
    if (!reader.Warning().empty()) std::cout << "TinyObjReader Warning: " << reader.Warning() << std::endl;

    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();
    bool enable_vertex_normals = (attrib.normals.size() > 0);

    auto object = make_shared<hittable_list>();

    // initialize num vertices & normals to have size 3
    std::vector<vec3> face_vertices(3);
    std::vector<vec3> face_normals(3);

    shared_ptr<material> __obj_mat;
    if (obj_material == nullptr) __obj_mat = make_shared<lambertian>(color::random()); //choose a random color
    else __obj_mat = obj_material;

    for (size_t s = 0; s < shapes.size(); s++) {
        
        // loop over the faces in the polygon
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            size_t fv = shapes[s].mesh.num_face_vertices[f];

            // loop over the verticies in the face
            for (size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                double vx = attrib.vertices[3*idx.vertex_index+0]*scale;
                double vy = attrib.vertices[3*idx.vertex_index+1]*scale;
                double vz = attrib.vertices[3*idx.vertex_index+2]*scale;
                face_vertices[v] = vec3(vx, vy, vz) + translation_vector;
                
                if (enable_vertex_normals) {
                    double nx = attrib.normals[3*idx.normal_index+0];
                    double ny = attrib.normals[3*idx.normal_index+1];
                    double nz = attrib.normals[3*idx.normal_index+2];
                    face_normals[v] = vec3(nx, ny, nz);
                }
            }
            if (enable_vertex_normals) {
                object->add(make_shared<triangle>(
                    face_vertices[0], face_vertices[1], face_vertices[2], 
                    face_normals[0], face_normals[1], face_normals[2],
                    __obj_mat));
            } else object->add(make_shared<triangle>(face_vertices[0], face_vertices[1], face_vertices[2], __obj_mat));
            index_offset += fv;
        }       
    }

    world.add(object);
    return true;
}

color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return color(0,0,0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}

hittable_list random_scene() {
    hittable_list world;

    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0,-1000,0), 1000, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

hittable_list load_mesh_scene() {
    hittable_list world;
    std::string obj_file = "obj_meshes/cow__1000_rotated.obj";
    
    // auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    // world.add(make_shared<sphere>(point3(0,-1001,0), 1000, ground_material));

    auto mesh_material = make_shared<lambertian>(color(0.25882, 0.5294, 0.96078));

    read_object_file(obj_file, 
        world, 
        5.0, 
        mesh_material, 
        vec3(0,0,10)
    );

    return world;
}

int main() {

    // Image

    const int image_width = 1280;
    const int image_height = 720;
    const double aspect_ratio = static_cast<double>(image_width) / static_cast<double>(image_height);
    const int samples_per_pixel = 10;
    const int max_depth = 50;

    // World

    auto world = load_mesh_scene();

    // Camera

    point3 lookfrom(0,0,0);
    point3 lookat(0,0,1);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.00001; //have large depth of field (i.e. no blur)

    camera cam(lookfrom, lookat, vup, 90, aspect_ratio, aperture, dist_to_focus);

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0,0,0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
