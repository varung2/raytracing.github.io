#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "rtweekend.h"

#include "hittable.h"

class triangle : public hittable {
	public:
		triangle() : enable_vertex_norms(false) {}
		triangle(const vec3& A, const vec3& B, const vec3& C, 
				 const vec3& A_norm, const vec3& B_norm, const vec3& C_norm,
				 const shared_ptr<material> m) : enable_vertex_norms(true) {
			corner = A;
			BA = B - A;
			CB = C - B;
			AC = A - C; 
			normal = unit_vector(cross(BA, CB));
			a = A; b = B; c = C;
			mat_ptr = m;

			a_norm = A_norm;
			b_norm = B_norm;
			c_norm = C_norm;

			inv_ABC_area = 1.0/dot(normal, cross(BA, -AC));
		}

		triangle(const vec3& A, const vec3& B, const vec3& C, const shared_ptr<material> m) : enable_vertex_norms(false) {//specify the points
			corner = A;
			BA = B - A;
			CB = C - B;
			AC = A - C; 
			normal = unit_vector(cross(BA, CB));
			a = A; b = B; c = C;
			mat_ptr = m;
		}

		virtual bool hit(const ray& r, double tmin, double tmax, hit_record &hit) const override;
		
	//member variables
	public:
		vec3 a, b, c; //points
		vec3 AC, BA, CB, normal, corner; //vectors
		shared_ptr<material> mat_ptr;
		// point normals
		vec3 a_norm, b_norm, c_norm;
	private:
		const bool enable_vertex_norms;
		double inv_ABC_area;
};

/** Function: hit(...)
	Description: Inherited hit function from parent. Checks if the ray has hit the object
	Inputs: 
		r 	  (type:ray) 	  	- casted ray that is checked against the geometry primitive
		t_min (type:double) 	- time minimum
		t_max (type:double) 	- time maximum
		rec   (type:hit_record) - record struct used to keep track of hit information if a hit is detected
	Outputs:
		bool  (type:boolean)	- returns if it hit something or not.
**/
bool triangle::hit(const ray& r, double tmin, double tmax, hit_record &rec) const {
	//calculate intersection using barycentric coordinates
	//First calculate the point of intersection on the plane
	double x = dot(normal, r.direction());
	if (fabs(x) < 0.001) return false;
	double t = (dot(corner - r.origin(), normal)/x);
	if (t < tmin || t > tmax) return false;
	vec3 p = r.at(t);

	vec3 AP = p - a;
	vec3 CP = p - c;
	vec3 BP = p - b;

	//topleft
	vec3 t0 = cross(BA, AP);
	double ABP_area = dot(normal, t0);
	if (ABP_area < 0) return false;
	//topright
	vec3 t1 = cross(CB, BP);
	double BCP_area = dot(normal, t1);
	if (BCP_area < 0) return false;
	//bottom
	vec3 t2 = cross(AC, CP);
	if (dot(normal, t2) < 0) return false;

	if (t < tmax && t > tmin) {
		rec.t = t;
		rec.p = p;
		rec.mat_ptr = mat_ptr;

		if (enable_vertex_norms) {
			vec3 smoothed_normal;
			double bary_a = BCP_area*inv_ABC_area;
			double bary_c = ABP_area*inv_ABC_area;
			smoothed_normal = unit_vector(bary_a*a_norm + (1.0 - bary_a - bary_c)*b_norm + bary_c*c_norm);
			// Bug fix 1: fixed color at the edges (normal can be facing away from the ray direction after smoothing) 
			//  - this causes a black color pixel because the scatter direction goes inside the object (only fixed lambertian shading)
			// Bug 2: getting dark artifacts on edges of metal objects
			//	- fix attempt 1: if the regular normal is facing the correct direction and the smooth is also, then use the smooth normal
			//					 otherwise, use the regular face normal (not working...)
			// bool front_face_regular = dot(r.direction(), normal) < 0;
			// bool front_face_smooth = dot(r.direction(), smoothed_normal) < 0;
			// smoothed_normal = (front_face_regular && front_face_smooth) ? smoothed_normal : normal;
			// smoothed_normal = (dot(r.direction(), smoothed_normal) < 0) ? smoothed_normal : normal;
			
			// if the smoothed normal is pointing in the right direction - use it
			if (dot(r.direction(), smoothed_normal) < 0) rec.normal = smoothed_normal;  
			else rec.set_face_normal(r, normal); //else set the normal to the regular face normal (w.r.t its orientation)
		} else rec.set_face_normal(r, normal);

		return true;
	}

	return false;
}

#endif