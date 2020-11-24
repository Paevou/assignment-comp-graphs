// Tampere University
// COMP.CE.430 Computer Graphics Coding Assignment 2020
//
// Write your name and student id here:
//   Kemppainen Paavo, H265380
//
// Mark here with an X which functionalities you implemented.
// Note that different functionalities are worth different amount of points.
//
// Name of the functionality      |Done| Notes
//-------------------------------------------------------------------------------
// example functionality          | X  | Example note: control this with var YYYY
// Mandatory functionalities ----------------------------------------------------
//   Perspective projection       |    | 
//   Phong shading                |    | 
//   Camera movement and rotation |    | 
//   Sharp shadows                |    | 
// Extra functionalities --------------------------------------------------------
//   Tone mapping                 |    | 
//   PBR shading                  |    | 
//   Soft shadows                 |    | 
//   Sharp reflections            |    | 
//   Glossy reflections           |    | 
//   Refractions                  |    | 
//   Caustics                     |    | 
//   SDF Ambient Occlusions       |    | 
//   Texturing                    |    | 
//   Simple game                  |    | 
//   Progressive path tracing     |    | 
//   Basic post-processing        |    | 
//   Advanced post-processing     |    | 
//   Screen space reflections     |    | 
//   Screen space AO              |    | 
//   Simple own SDF               |    | 
//   Advanced own SDF             |    | 
//   Animated SDF                 |    | 
//   Other?                       |    | 


#ifdef GL_ES
precision mediump float;
#endif


// Defines for functionalities that cannot be on at the same time
#define SHARP_SHADOW false
#define SOFT_SHADOW false

#define SHARP_REFLECTION false
#define GLOSSY_REFLECTION false

#define REFRACTIONS true



#define PI 3.14159265359
#define EPSILON 0.00001

// These definitions are tweakable.

/* Minimum distance a ray must travel. Raising this value yields some performance
 * benefits for secondary rays at the cost of weird artefacts around object
 * edges.
 */
#define MIN_DIST 0.08
/* Maximum distance a ray can travel. Changing it has little to no performance
 * benefit for indoor scenes, but useful when there is nothing for the ray
 * to intersect with (such as the sky in outdoors scenes).
 */
#define MAX_DIST 20.0
/* Maximum number of steps the ray can march. High values make the image more
 * correct around object edges at the cost of performance, lower values cause
 * weird black hole-ish bending artefacts but is faster.
 */
#define MARCH_MAX_STEPS 128
/* Typically, this doesn't have to be changed. Lower values cause worse
 * performance, but make the tracing stabler around slightly incorrect distance
 * functions.
 * The current value merely helps with rounding errors.
 */
#define STEP_RATIO 0.999
/* Determines what distance is considered close enough to count as an
 * intersection. Lower values are more correct but require more steps to reach
 * the surface
 */
#define HIT_RATIO 0.001

// Resolution of the screen
uniform vec2 u_resolution;

// Mouse coordinates
uniform vec2 u_mouse;

// Time since startup, in seconds
uniform float u_time;

// Values are approximated with exercise 5 values
// Ambient constant
const float AMBIENT_STRENGTH = 0.6;
// Diffuse constant
const vec3 DIFFUSE = vec3(0.7, 0.7, 0.3);
// Specular constant
const vec3 SPECULAR = vec3(0.7, 0.6, 0.7);
const float DIFFUSE_INTENSITY = 0.1;
const float SPECULAR_INTENSITY = 0.2;
const float SHININESS = 0.2;

const float light_bulb_radius = 0.5;
const vec3 lamp_pos = vec3(0.0, 2.1, 3.0);
const float air_refraction_term = 1.0;

float displace(vec3 p);

float sharp_shadow(vec3 object_point, vec3 normal, vec3 light_dir, vec3 light_point);
float soft_shadow(vec3 object_point,  vec3 normal, vec3 light_dir, vec3 light_point);
vec3 getConeSampleShadow(vec3 direction, float coneAngle, vec3 perpendicular, in vec2 in_seed, out vec2 out_seed, float dist);
float rand(in vec2 in_seed, out vec2 out_seed);
mat4 rotationMatrix(vec3 axis, float angle);
vec3 getConeSampleReflection(vec3 direction, float coneAngle, vec3 perpendicular, in vec2 in_seed, out vec2 out_seed, float dist);


vec3 render_reflect(vec3 o, vec3 v);
vec3 render_refraction(vec3 o, vec3 v);
vec3 sharp_reflection(vec3 object_point, vec3 reflect_dir);

struct material
{
    // The color of the surface
    vec4 color;
    // You can add your own material features here!
    float reflection_term;
    float refraction_term;
    float transparency_term;
};
vec3 glossy_reflection(vec3 object_point, vec3 reflect_dir, material o_mat);
vec3 refraction(vec3 o, vec3 v, vec3 n, material o_mat);
vec3 refraction_vector(vec3 n, vec3 v, float n1, float n2);

// Good resource for finding more building blocks for distance functions:
// https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

/* Basic box distance field.
 *
 * Parameters:
 *  p   Point for which to evaluate the distance field
 *  b   "Radius" of the box
 *
 * Returns:
 *  Distance to the box from point p.
 */
float box(vec3 p, vec3 b)
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

// Simple cylinder distance field
float cylinder(vec3 p, vec3 r) {
    // p is the place where the object is going to be and r is the dimensions
    vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r.x,r.y);
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float fractal(vec3 p, vec3 r) {
    float boundary = 50;
    float len = 0;
    vec2 z = vec2(0.0);
    for(int i = 0; i < 100; i++ ) {
        z = vec2(z.x * z.x - z.y * z.y, 2.0 * z.x * z.y ) + p.xy;
        if( dot(z,z) > (boundary * boundary)) break;
        len += 1.0;
    }

    if( len > 99 ) {
        return 0.0;
    }
    float s_len = len - log2(log2(dot(z,z))) + 4.0;

    float a_len = smoothstep( -0.1, 0.0, sin(0.5 * 6.2831 * u_time));
    return len;

    vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r.x,r.y);
    float d2 = displace(p);
    return min(max(d.x + d2, d.y + d2),d.x);
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float displace(vec3 p) {
    return sin(20.0*p.x)*sin(20.0*p.y)*sin(20.0*p.z);
}

/* Rotates point around origin along the X axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_x(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        p.x,
        c*p.y-s*p.z,
        s*p.y+c*p.z
    );
}

/* Rotates point around origin along the Y axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_y(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        c*p.x+s*p.z,
        p.y,
        -s*p.x+c*p.z
    );
}

/* Rotates point around origin along the Z axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_z(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        c*p.x-s*p.y,
        s*p.x+c*p.y,
        p.z
    );
}

/* Each object has a distance function and a material function. The distance
 * function evaluates the distance field of the object at a given point, and
 * the material function determines the surface material at a point.
 */

float blob_distance(vec3 p)
{
    vec3 q = p - vec3(-0.5, -2.2 + abs(sin(u_time*3.0)), 2.0);
    return length(q) - 0.8 + sin(10.0*q.x)*sin(10.0*q.y)*cos(10.0*q.z)*0.07;
}

material blob_material(vec3 p)
{
    material mat;
    mat.reflection_term = 0.0;
    mat.refraction_term = 0.0;
    mat.transparency_term = 0.5;
    mat.color = vec4(1.0, 0.5, 0.3, 0.0);
    return mat;
}

float sphere_distance(vec3 p)
{
    return length(p - vec3(1.5, -1.8, 2.0)) - 1.2;
}

material sphere_material(vec3 p)
{
    material mat;
    mat.reflection_term = 0.2;
    mat.refraction_term = 1.3;
    mat.transparency_term = 1.0;
    mat.color = vec4(0.1725, 0.4392, 0.1882, 1.0);
    return mat;
}

float room_distance(vec3 p)
{
    return max(
        -box(p-vec3(0.0,3.0,3.0), vec3(0.5, 0.5, 0.5)),
        -box(p-vec3(0.0,0.0,0.0), vec3(3.0, 3.0, 6.0))
    );
}

material room_material(vec3 p)
{
    material mat;
    mat.reflection_term = 0.1;
    mat.refraction_term = 0.0;
    mat.transparency_term = 0.0;
    mat.color = vec4(1.0, 1.0, 1.0, 1.0);
    if(p.x <= -2.98) mat.color.rgb = vec3(1.0, 0.0, 0.0);
    else if(p.x >= 2.98) mat.color.rgb = vec3(0.0, 1.0, 0.0);
    return mat;
}

float crate_distance(vec3 p)
{
    return box(rot_y(p-vec3(-1,-1,5), u_time), vec3(1, 2, 1));
}

material crate_material(vec3 p)
{
    material mat;
    mat.color = vec4(1.0, 1.0, 1.0, 1.0);
    mat.reflection_term = 0.7;
    mat.refraction_term = 1.1;
    mat.transparency_term = 0.0;
    vec3 q = rot_y(p-vec3(-1,-1,5), u_time) * 0.98;
    if(fract(q.x + floor(q.y*2.0) * 0.5 + floor(q.z*2.0) * 0.5) < 0.5)
    {
        mat.color.rgb = vec3(0.0, 1.0, 1.0);
    }
    return mat;
}

float light_distance(vec3 p)
{
    return length(p - lamp_pos) - light_bulb_radius;
}

material light_material(vec3 p)
{
    material mat;
    mat.reflection_term = 0.0;
    mat.refraction_term = 0.0;
    mat.transparency_term = 0.0;
    mat.color = vec4(0.851, 1.0, 0.0, 1.0);
    return mat;
}

float cylinder_distance(vec3 p) {
    return cylinder(p - vec3(1.5, -1.8, 5.0), vec3(1.0, 1.5, 1.0));
}

material cylinder_material(vec3 p ) {
    material mat;
    mat.reflection_term = 0.2;
    mat.refraction_term = 1.3;
    mat.transparency_term = 0.0;
    mat.color = vec4(0.4, 0.8275, 0.4275, 1.0);
    return mat;
}

float fractal_distance(vec3 p) {
    return fractal(p - vec3(-1.0, -2.0, 3.0), vec3(1.0, 1.5, 1.0));
}

material fractal_material(vec3 p ) {
    material mat;
    mat.reflection_term = 0.2;
    mat.refraction_term = 1.3;
    mat.transparency_term = 0.0;
    mat.color = vec4(0.4, 0.8275, 0.4275, 1.0);
    return mat;
}


/* The distance function collecting all others.
 *
 * Parameters:
 *  p   The point for which to find the nearest surface
 *  mat The material of the nearest surface
 *
 * Returns:
 *  The distance to the nearest surface.
 */
float map(
    in vec3 p,
    out material mat
){
    float min_dist = MAX_DIST*2.0;
    float dist = 0.0;

    // dist = blob_distance(p);
    // if(dist < min_dist) {
    //     mat = blob_material(p);
    //     min_dist = dist;
    // }

    dist = room_distance(p);
    if(dist < min_dist) {
        mat = room_material(p);
        min_dist = dist;
    }

    dist = crate_distance(p);
    if(dist < min_dist) {
        mat = crate_material(p);
        min_dist = dist;
    }

    dist = sphere_distance(p);
    if(dist < min_dist) {
        mat = sphere_material(p);
        min_dist = dist;
    }

    // Add your own objects here!

    dist = light_distance(p);
    if(dist < min_dist) {
        mat = light_material(p);
        min_dist = dist;
    }

    dist = cylinder_distance(p);
    if(dist < min_dist) {
        mat = cylinder_material(p);
        min_dist = dist;
    }

    dist = fractal_distance(p);
    if(dist < min_dist) {
        mat = fractal_material(p);
        min_dist = dist;
    }
    
    return min_dist;
}

/* Calculates the normal of the surface closest to point p.
 *
 * Parameters:
 *  p   The point where the normal should be calculated
 *  mat The material information, produced as a byproduct
 *
 * Returns:
 *  The normal of the surface.
 *
 * See https://www.iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
 * if you're interested in how this works.
 */
vec3 normal(vec3 p, out material mat)
{
    const vec2 k = vec2(1.0, -1.0);
    return normalize(
        k.xyy * map(p + k.xyy * EPSILON, mat) +
        k.yyx * map(p + k.yyx * EPSILON, mat) +
        k.yxy * map(p + k.yxy * EPSILON, mat) +
        k.xxx * map(p + k.xxx * EPSILON, mat)
    );
}

/* Finds the closest intersection of the ray with the scene.
 *
 * Parameters:
 *  o           Origin of the ray
 *  v           Direction of the ray
 *  max_dist    Maximum distance the ray can travel. Usually MAX_DIST.
 *  p           Location of the intersection
 *  n           Normal of the surface at the intersection point
 *  mat         Material of the intersected surface
 *  inside      Whether we are marching inside an object or not. Useful for
 *              refractions.
 *
 * Returns:
 *  true if a surface was hit, false otherwise.
 */
bool intersect(
    in vec3 o,
    in vec3 v,
    in float max_dist,
    out vec3 p,
    out vec3 n,
    out material mat,
    bool inside
) {
    float t = MIN_DIST;
    float dir = inside ? -1.0 : 1.0;
    bool hit = false;

    for(int i = 0; i < MARCH_MAX_STEPS; ++i)
    {
        p = o + t * v;
        float dist = dir * map(p, mat);
        
        hit = abs(dist) < HIT_RATIO * t;

        if(hit || t > max_dist) break;

        t += dist * STEP_RATIO;
    }

    n = normal(p, mat);

    return hit;
}

/* Calculates the color of the pixel, based on view ray origin and direction.
 *
 * Parameters:
 *  o   Origin of the view ray
 *  v   Direction of the view ray
 *
 * Returns:
 *  Color of the pixel.
 */
vec3 render(vec3 o, vec3 v)
{
    // This lamp is positioned at the hole in the roof.

    vec3 p, n;
    material mat;

    // Compute intersection point along the view ray.
    intersect(o, v, MAX_DIST, p, n, mat, false);

    // Add some lighting code here! 

    // Phong lighting   
    vec3 light_dir = normalize(lamp_pos - p);
    v = normalize(v);
    vec3 diffuse = DIFFUSE * DIFFUSE_INTENSITY * max(dot(light_dir, n), 0.0);
    vec3 reflect_ld = 2.0 * dot(-light_dir, n) * n + light_dir;
    vec3 specular = SPECULAR_INTENSITY * SPECULAR * pow(max(dot(reflect_ld, v), 0.0), SHININESS);

    mat.color.rgb = AMBIENT_STRENGTH * mat.color.rgb + diffuse + specular;
    // Come out of the surface
    p += 0.01*n;
    if( SHARP_SHADOW ) {
        // Sharp shadow        
        mat.color.rgb *= sharp_shadow(p, n, light_dir, lamp_pos);
    } else if (SOFT_SHADOW) {
        // Soft shadow  
        mat.color.rgb *= soft_shadow(p, n, light_dir, lamp_pos);           
    }
    
    vec3 reflect_dir = reflect(v, n);
    vec3 reflect_color = vec3(0.0);
    if(SHARP_REFLECTION ) {
        reflect_color = sharp_reflection(p, reflect_dir);
    } else if (GLOSSY_REFLECTION) {
        reflect_color = glossy_reflection(p, reflect_dir, mat);
    }
    vec3 refract_color = vec3(0.0);
    p -= 0.01*n;
    if(REFRACTIONS) {
        refract_color = refraction(p, v, n, mat);
    }
    
    
    mat.color.rgb += reflect_color * mat.reflection_term;
    if(mat.transparency_term != 0.0 && refract_color != vec3(0.0)) {
        mat.color.rgb *= 1.0 - mat.transparency_term;
        mat.color.rgb += refract_color * mat.transparency_term;
    }

    return mat.color.rgb;
}

float sharp_shadow(vec3 object_point, vec3 normal, vec3 light_dir, vec3 light_point) {
    vec3 p, n;
    material mat;

    float dist = sqrt( pow(object_point.x - light_point.x,2.0) 
        + pow(object_point.y - light_point.y,2.0) 
        + pow(object_point.z - light_point.z,2.0));
    object_point += 0.01*normal;
    bool intersect = intersect(object_point, light_dir, dist, p, n, mat, false);
    
    //if( p != light_point && intersect) {
    if( mat != light_material(vec3(0.0))) {
        return 0.5;
    } else {
        return 1.0;
    }    
}

float soft_shadow(vec3 object_point, vec3 normal, vec3 light_dir, vec3 light_point) {
    // Source which was used as a base for soft shadows
    //https://medium.com/@alexander.wester/ray-tracing-soft-shadows-in-real-time-a53b836d123b
    
    vec3 perpendicular = cross(light_dir, vec3(0.0,1.0,0.0));
    if(perpendicular.x == 0.0
        && perpendicular.y == 0.0
        && perpendicular.z == 0.0) {
        perpendicular = vec3(1.0,0.0,0.0);
    } 

    // Vector to the edge of the light bulb
    vec3 vec_light_edge = normalize((light_point + perpendicular * light_bulb_radius) - object_point);
    // * 2 so we get both halves of the circle plane
    float cone_angle = acos(dot(light_dir, vec_light_edge)) * 2.0;

    const int num_of_rays = 2;
    int num_of_hits = 0;
    object_point += 0.01*normal;    
    vec2 in_seed = gl_FragCoord.xy/u_resolution.xy;
    vec2 out_seed;
    for(int i = 0; i < num_of_rays; i++ ) {
        vec3 p, n;
        material mat;
        float dist = sqrt( pow(object_point.x - light_point.x,2.0) 
        + pow(object_point.y - light_point.y,2.0) 
        + pow(object_point.z - light_point.z,2.0));

        //vec3 t = getConeSample(seed, -light_dir, cone_angle); 
        vec3 t = getConeSampleShadow(light_dir, cone_angle, perpendicular, in_seed, out_seed, dist);
        in_seed = out_seed;
        vec3 ray_dir = light_dir * dist + t;
        ray_dir = t;
        //ray_dir = vec3(1.0,0.0,0.0);

        bool hit = intersect(object_point, ray_dir, dist, p, n, mat, false);
        if(mat == light_material(vec3(0.0))) {
            num_of_hits += 1;
        }        
    }    
    
    return float(num_of_hits)/float(num_of_rays);
}

vec3 getConeSampleShadow(vec3 direction, float coneAngle, vec3 perpendicular, in vec2 in_seed, out vec2 out_seed, float dist) {
    
    float r = light_bulb_radius * sin(rand(in_seed, out_seed));
    if(r == 0.0 ) {
        return direction;
    }
    float angle_around_axis = acos(rand(in_seed, out_seed));
    vec4 v = normalize(vec4(perpendicular, 1.0));
    mat4 rot_mat = rotationMatrix(direction, angle_around_axis);
    v = rot_mat * v * r;
    vec3 ret = dist * direction + v.xyz;
    ret = normalize(ret);
    return ret;
}

vec3 getConeSampleReflection(vec3 direction, float coneAngle, vec3 perpendicular, in vec2 in_seed, out vec2 out_seed, float dist) {
    float r = 1.0 * sin(rand(in_seed, out_seed)) * sin(coneAngle);
    if(r == 0.0 ) {
        return direction;
    }
    float angle_around_axis = acos(rand(in_seed, out_seed));
    vec4 v = normalize(vec4(perpendicular, 1.0));
    mat4 rot_mat = rotationMatrix(direction, angle_around_axis);
    v = rot_mat * v * r;
    vec3 ret = dist * direction + v.xyz;
    ret = normalize(ret);
    return ret;
}

float rand(in vec2 in_seed, out vec2 out_seed){
    // 'Random' function 
    out_seed = in_seed * 2.0;
    return fract(sin(dot(in_seed.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}

mat4 rotationMatrix(vec3 axis, float angle)
{
    // Source: http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

vec3 render_reflect(vec3 o, vec3 v)
{
    vec3 p, n;
    material mat;
    intersect(o, v, MAX_DIST, p, n, mat, false);  
    vec3 light_dir = normalize(lamp_pos - p);
    v = normalize(v);
    vec3 diffuse = DIFFUSE * DIFFUSE_INTENSITY * max(dot(light_dir, n), 0.0);
    vec3 reflect_ld = 2.0 * dot(-light_dir, n) * n + light_dir;
    vec3 specular = SPECULAR_INTENSITY * SPECULAR * pow(max(dot(reflect_ld, v), 0.0), SHININESS);

    mat.color.rgb = AMBIENT_STRENGTH * mat.color.rgb + diffuse + specular;
    // Come out of the surface
    p += 0.01*n;
    if( SHARP_SHADOW ) {
        // Sharp shadow        
        float sharp_shadow_term = sharp_shadow(p, n, light_dir, lamp_pos);
        mat.color.rgb *= sharp_shadow_term;
    } else if (SOFT_SHADOW) {
        // Soft shadow 
        float sharp_shadow_term = sharp_shadow(p, n, light_dir, lamp_pos);
        if( sharp_shadow_term == 1.0) {
            mat.color.rgb *= sharp_shadow_term;
        } else {
            mat.color.rgb *= soft_shadow(p, n, light_dir, lamp_pos); 
        }   
        
    }
    return mat.color.rgb;
}

vec3 render_refraction(vec3 o, vec3 v)
{
    // This lamp is positioned at the hole in the roof.

    vec3 p, n;
    material mat;

    // Compute intersection point along the view ray.
    intersect(o, v, MAX_DIST, p, n, mat, false);

    // Add some lighting code here! 

    // Phong lighting   
    vec3 light_dir = normalize(lamp_pos - p);
    v = normalize(v);
    vec3 diffuse = DIFFUSE * DIFFUSE_INTENSITY * max(dot(light_dir, n), 0.0);
    vec3 reflect_ld = 2.0 * dot(-light_dir, n) * n + light_dir;
    vec3 specular = SPECULAR_INTENSITY * SPECULAR * pow(max(dot(reflect_ld, v), 0.0), SHININESS);

    mat.color.rgb = AMBIENT_STRENGTH * mat.color.rgb + diffuse + specular;
    // Come out of the surface
    p += 0.01*n;
    if( SHARP_SHADOW ) {
        // Sharp shadow        
        mat.color.rgb *= sharp_shadow(p, n, light_dir, lamp_pos);
    } else if (SOFT_SHADOW) {
        // Soft shadow  
        mat.color.rgb *= soft_shadow(p, n, light_dir, lamp_pos);           
    }
    
    vec3 reflect_dir = reflect(v, n);
    vec3 reflect_color = vec3(0.0);
    if(SHARP_REFLECTION ) {
        reflect_color = sharp_reflection(p, reflect_dir);
    } else if (GLOSSY_REFLECTION) {
        reflect_color = glossy_reflection(p, reflect_dir, mat);
    }

    mat.color.rgb += reflect_color * mat.reflection_term;

    return mat.color.rgb;
}

vec3 sharp_reflection(vec3 object_point, vec3 reflect_dir) {
    vec3 p, n;
    material mat;
    bool hit = intersect(object_point, reflect_dir, MAX_DIST, p, n, mat, false );
    if(hit) {
        return render_reflect(object_point, reflect_dir);
    }
    return vec3(0.0);
}

vec3 glossy_reflection(vec3 object_point, vec3 reflect_dir, material o_mat) {

    vec3 perpendicular = cross(reflect_dir, vec3(0.0,1.0,0.0));

    if(perpendicular.x == 0.0
        && perpendicular.y == 0.0
        && perpendicular.z == 0.0) {
        perpendicular = vec3(1.0,0.0,0.0);
    }
    vec2 in_seed = gl_FragCoord.xy/u_resolution.xy;
    vec2 out_seed;

    vec3 p, n;
    material mat;
    const int num_of_rays = 2;
    int num_of_hits = 0;
    float angle = PI/16.0 * o_mat.reflection_term;
    vec3 color = vec3(0.0);
    for(int i = 0 ; i < num_of_rays; i++) {
        reflect_dir = getConeSampleReflection(reflect_dir,angle, perpendicular, in_seed, out_seed, 1.0);
        bool hit = intersect(object_point, reflect_dir, MAX_DIST, p, n, mat, false );
        if(hit) {
            num_of_hits += 1;
            color += render_reflect(object_point, reflect_dir)*(mat.reflection_term+1.5);
        }
    } 
    if(num_of_hits == 0 ) {
        return vec3(1.0);
    }
    return color/float(num_of_hits);
}

vec3 refraction(vec3 o, vec3 v, vec3 n, material o_mat) {
    // Some source: https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
    if(o_mat.refraction_term == 0.0) {
        return vec3(0.0);
    }
    //vec3 ref_v = refraction_vector(n, v, air_refraction_term,
    //o_mat.refraction_term);
    vec3 ref_v = refract(v, n, air_refraction_term/o_mat.refraction_term);
    vec3 p;
    o += 0.01 * -n;
    material mat;
    intersect(o, ref_v, MAX_DIST, p, n, mat, true);

    //ref_v = refraction_vector(n, v, o_mat.refraction_term, air_refraction_term);
    ref_v = refract(v, n, o_mat.refraction_term/air_refraction_term);
    p += 0.01*n;
    vec3 color = render_refraction(p, -ref_v);
    return color;
}

vec3 refraction_vector(vec3 n, vec3 v, float n1, float n2) {
    float n12 = n1/n2;
    float cosV = -dot(n, v);
    float sinT2 = n12 * n12 * (1.0 - cosV * cosV);
    if(sinT2 > 1.0) {
        //TODO: Some illegal value should be returned
        //return vec3(0.0);
    }
    float cosT = sqrt(1.0 - sinT2);
    return n12 * v + (n12 * cosV - cosT ) * n;
}

void main()
{    
    // This is the position of the pixel in normalized device coordinates.
    vec2 uv = (gl_FragCoord.xy/u_resolution)*2.0-1.0;
    // Calculate aspect ratio
    float aspect = u_resolution.x/u_resolution.y;

    // Modify these two to create perspective projection!
    //TODO: Why do parts of the blod dissappear
    // Origin of the view ray       
    //vec3 o = vec3(2.96*vec2(uv.x * aspect, uv.y), -2.0);
    vec3 o = vec3(0,0,-0.8);
    // Direction of the view ray
    //vec3 v = vec3(0,0,1);
    vec3 v = vec3(uv.x, uv.y, 0) - o;

    // Camera movement
    // vec3 camera_pos = o;
    // TODO: Possibly have to change to +o, as then camera points towards
    // positive z-axis
    // vec3 camera_dir = normalize(-o);
    // vec3 camera_up = normalize(cross(o, vec3(0.0, 1.0, -1.0)));
    // vec3 camera_right = normalize(cross(camera_dir, camera_up));
    // camera_up = cross(camera_right, camera_dir);
    
    // TODO: Ask about wheter camera movement is in fragment shader or the
    // vertex shader
    // float norm_mouse_x = -1.0 + 2.0 * u_mouse.x / u_resolution.x;
    // float norm_mouse_y = -1.0 + 2.0 * u_mouse.y / u_resolution.y;
    // v = rot_y(v, norm_mouse_x);
    // v = rot_x(v, -norm_mouse_y);
    // o = vec3(1.0*sin(u_time), 1.0*cos(u_time),1.0*sin(u_time));
    // v = rot_z(v, u_mouse.y/u_resolution.x);
    
    gl_FragColor = vec4(render(o, v), 1.0);
}
