#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "udray.h"
#include "glm.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

extern Camera *ray_cam;       // camera info
extern int image_i, image_j;  // current pixel being shaded
extern bool wrote_image;      // has the last pixel been shaded?

// reflection/refraction recursion control

extern int maxlevel;          // maximum depth of ray recursion 
extern double minweight;      // minimum fractional contribution to color
extern double rayeps;

// these describe the scene

extern vector < GLMmodel * > model_list;
extern vector < Sphere * > sphere_list;
extern vector < Light * > light_list;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// intersect a ray with the entire scene (.obj models + spheres)

// x, y are in pixel coordinates with (0, 0) the upper-left hand corner of the image.
// color variable is result of this function--it carries back info on how to draw the pixel

void trace_ray(int level, double weight, Ray *ray, Vect color)
{
	Intersection *nearest_inter = NULL;
	Intersection *inter = NULL;
	int i;

	// test for intersection with all .obj models

	for (i = 0; i < model_list.size(); i++) {
		inter = intersect_ray_glm_object(ray, model_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	// test for intersection with all spheres

	for (i = 0; i < sphere_list.size(); i++) {
		inter = intersect_ray_sphere(ray, sphere_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	// "color" the ray according to intersecting surface properties

	// choose one of the simpler options below to debug or preview your scene more quickly.
	// another way to render faster is to decrease the image size.

	if (nearest_inter) {
		//shade_ray_false_color_normal(nearest_inter, color);
		//    shade_ray_intersection_mask(color);  
	//	shade_ray_diffuse(ray, nearest_inter, color);
	//	shade_ray_local(ray, nearest_inter, color);
		shade_ray_recursive(level, weight, ray, nearest_inter, color);
	}

	// color the ray using a default

	else
		shade_ray_background(ray, color); 
}

//----------------------------------------------------------------------------

// test for ray-sphere intersection; return details of intersection if true

Intersection *intersect_ray_sphere(Ray *ray, Sphere *S)
{
	// FILL IN CODE (line below says "no" for all spheres, so replace it)
	Vect    a,I;        // triangle vectors
	double   b,c,d;           // params to calc ray-plane intersect
	double t,deta,t1,t2;
	Intersection *inter;
	//VectPrint(ray->dir);

	VectSub(ray->orig, S->P, a); //vector "-"
	c=2.0*VectDotProd(a, ray->dir);
	//VectPrint(ray->orig);
	//VectPrint(S->P);
	//VectPrint(a);
	//VectPrint(ray->dir);
	//cout << "c="<< c << endl;
	b=VectDotProd(ray->dir, ray->dir);
	d = VectDotProd(a, a) - S->radius * S->radius;
	//cout << "r=" << S->radius << endl;
	deta = c*c - 4*b*d;
	//cout << "deta=" << deta << endl;
	inter = make_intersection();
	if (deta < 0)
		return NULL;
	if (deta >= 0)
	{
		t1 = (-c + sqrt(deta)) / (2 * b);
		t2= (-c - sqrt(deta)) / (2 * b);
		//cout << "t1=" << t1 << " " << "t2=" << t2 << endl;
		if (t1 > t2)
			if (t2 >= 0)
			inter->t = t2;
		else
			if (t1 >= 0)
			inter->t = t1;
	}
	//cout << inter->t << endl;
	VectAddS(inter->t, ray->dir, ray->orig, I);
	VectCopy(inter->P, I);
	return inter;
	// return NULL;
}
//----------------------------------------------------------------------------

// only local, ambient + diffuse lighting (no specular, shadows, reflections, or refractions)

void shade_ray_diffuse(Ray *ray, Intersection *inter, Vect color)
{
	Vect L,R_prim;
	Vect normal;
	Vect a,c,V,V_unit;
	double phong;
	double diff_factor=1;
	double b;
	double diffuse_R,diffuse_G,diffuse_B;
	double specular_R,specular_G,specular_B;

	// iterate over lights

	for (int i = 0; i < light_list.size(); i++) {

		// AMBIENT

		color[R] += inter->surf->amb[R] * light_list[i]->amb[R];
		color[G] += inter->surf->amb[G] * light_list[i]->amb[G];
		color[B] += inter->surf->amb[B] * light_list[i]->amb[B];
		//cout<<color[G]<<endl;
		// DIFFUSE
		VectSub(ray->orig,inter->P,V);
		VectNormalize(V, V_unit);
		VectSub(light_list[i]->P, inter->P, a);
		VectNormalize(a, L);
		VectNormalize(inter->N, normal);
		b=VectDotProd(normal, L);
		VectNumber(2*b,normal,c);
		VectSub(c,L,R_prim);
		phong=VectDotProd(V_unit, R_prim);
		//cout<<phong<<endl;
		if (phong<=0)
		{
			specular_R=specular_G=specular_B=0;
		}
		else
		{
		//for (int j=0; j < inter->surf->spec_exp;j++)
		//{
		//		phong*=phong;
		//		cout<<phong<<endl;
		//}
		phong=pow(phong,inter->surf->spec_exp);
		//cout<<"inter->surf->spec[R]="<<inter->surf->spec[R]<<endl;
		//cout<<"phong="<<phong<<endl;
		//cout<<"light_list[i]->spec[R]="<<light_list[i]->spec[R]<<endl;
		//cout<<"light_list[i]->amb[R]="<<light_list[i]->amb[R]<<endl;
		//cout<<"light_list[i]->diff[R]="<<light_list[i]->diff[R]<<endl;
		//specular_R=inter->surf->spec[R]*phong*light_list[i]->spec[R]*diff_factor*1.2;
		//specular_G=inter->surf->spec[G]*phong*light_list[i]->spec[G]*diff_factor*1.2;
		//specular_B=inter->surf->spec[B]*phong*light_list[i]->spec[B]*diff_factor*1.2;

		//specular_R=inter->surf->spec[R]*phong*diff_factor;
		//specular_G=inter->surf->spec[G]*phong*diff_factor;
		//specular_B=inter->surf->spec[B]*phong*diff_factor;

		}

		if (b<=0)
		{
			diffuse_R=diffuse_G=diffuse_B=0;
		}
		else
		{
			diffuse_R=inter->surf->diff[R]*b*light_list[i]->diff[R]*diff_factor*1.0;
			diffuse_G=inter->surf->diff[G]*b*light_list[i]->diff[G]*diff_factor*1.0;
			diffuse_B=inter->surf->diff[B]*b*light_list[i]->diff[B]*diff_factor*1.0;
		}
		color[R] = diffuse_R+color[R]*(1-diff_factor);
		color[G] = diffuse_G+color[G]*(1-diff_factor);
		color[B] = diffuse_B+color[B]*(1-diff_factor);
		//cout<<color[G]<<endl;
		// FILL IN CODE

	}

	// clamp color to [0, 1]

	VectClamp(color, 0, 1);
}

//----------------------------------------------------------------------------

// same as shade_ray_diffuse(), but add specular lighting + shadow rays (i.e., full Phong illumination model)

void shade_ray_local(Ray *ray, Intersection *inter, Vect color)
{
	
	Vect L,R_prim;
	Vect normal;
	Vect a,c,V,V_unit;
	double phong;
	double specular_factor=0.9;
	double b;
	double diffuse_R,diffuse_G,diffuse_B;
	double specular_R,specular_G,specular_B;

	// iterate over lights

	for (int i = 0; i < light_list.size(); i++) {

		// AMBIENT

		color[R] += inter->surf->amb[R] * light_list[i]->amb[R];
		color[G] += inter->surf->amb[G] * light_list[i]->amb[G];
		color[B] += inter->surf->amb[B] * light_list[i]->amb[B];
		//cout<<color[G]<<endl;
		// DIFFUSE
		VectSub(ray->orig,inter->P,V);
		VectNormalize(V, V_unit);
		VectSub(light_list[i]->P, inter->P, a);
		VectNormalize(a, L);
		VectNormalize(inter->N, normal);
		b=VectDotProd(normal, L);
		VectNumber(2*b,normal,c);
		VectSub(c,L,R_prim);
		phong=VectDotProd(V_unit, R_prim);
		//cout<<phong<<endl;
		if (phong<=0)
		{
			specular_R=specular_G=specular_B=0;
		}
		else
		{
		//for (int j=0; j < inter->surf->spec_exp;j++)
		//{
		//		phong*=phong;
		//		cout<<phong<<endl;
		//}
		phong=pow(phong,inter->surf->spec_exp);
		//cout<<"inter->surf->spec[R]="<<inter->surf->spec[R]<<endl;
		//cout<<"phong="<<phong<<endl;
		//cout<<"light_list[i]->spec[R]="<<light_list[i]->spec[R]<<endl;
		//cout<<"light_list[i]->amb[R]="<<light_list[i]->amb[R]<<endl;
		//cout<<"light_list[i]->diff[R]="<<light_list[i]->diff[R]<<endl;
		//specular_R=inter->surf->spec[R]*phong*light_list[i]->spec[R]*diff_factor*1.2;
		//specular_G=inter->surf->spec[G]*phong*light_list[i]->spec[G]*diff_factor*1.2;
		//specular_B=inter->surf->spec[B]*phong*light_list[i]->spec[B]*diff_factor*1.2;
		specular_R=inter->surf->spec[R]*phong*specular_factor;
		specular_G=inter->surf->spec[G]*phong*specular_factor;
		specular_B=inter->surf->spec[B]*phong*specular_factor;

		}

		color[R] += specular_R;
		color[G] += specular_G;
		color[B] += specular_B;
		//cout<<color[G]<<endl;
		// FILL IN CODE

	}

	// clamp color to [0, 1]

	VectClamp(color, 0, 1);
	// FILL IN CODE 
}

//----------------------------------------------------------------------------

// full shading model: ambient/diffuse/specular lighting, shadow rays, recursion for reflection, refraction

// level = recursion level (only used for reflection/refraction)

void shade_ray_recursive(int level, double weight, Ray *ray, Intersection *inter, Vect color)
{
	Surface *surf;
	int i;
	Vect Modify_P;
	Vect N_ray_dir;
	Vect L;
	Vect Zero_vect;
	Vect L_vert;
	Vect L_hori;
	Vect Reflect_vect;
	Vect Refract_vect;
	Ray*  inter_ray=make_ray();
	double Reflect_factor;
	double Refract_factor;
	double theta1;
	double theta2;
	double t=0;


	// initialize color to Phong reflectance model

	//shade_ray_local(ray, inter, color);
		shade_ray_diffuse(ray, inter, color);
		shade_ray_local(ray, inter, color);

	// if not too deep, recurse
	surf=inter->surf;
    VectClamp(Zero_vect,0.0,0.0);
//cout<<level<<" "<<maxlevel<<endl;
	if (level + 1 < maxlevel) {

		// add reflection component to color

		if (surf->reflectivity * weight > minweight) {

			// FILL IN CODE
			    VectAddS(rayeps,inter->N,inter->P,Modify_P);

			//reflect dir
			VectNegate(ray->dir,N_ray_dir);
			Reflect_factor=VectDotProd(inter->N,N_ray_dir);
			VectAddS(Reflect_factor,inter->N,Zero_vect,L_vert);
			VectSub(N_ray_dir,L_vert,L_hori);
			VectSub(L_vert,L_hori,Reflect_vect);

			//copy
			VectCopy(inter_ray->orig,Modify_P);
			VectCopy(inter_ray->dir,Reflect_vect);


			trace_ray(level+1,weight,inter_ray,color);

		}

		// add refraction component to color

		if (surf->transparency * weight > minweight) {

			// GRAD STUDENTS -- FILL IN CODE
					        //add epx
		        VectAddS(-rayeps,inter->N,inter->P,Modify_P);

				VectCopy(L,ray->dir);

			    VectNegate(ray->dir,N_ray_dir);

				//cos(\theta_{1})
			    Refract_factor=VectDotProd(inter->N,N_ray_dir);
				//theta_{1}
				theta1=acos(Refract_factor);
				theta2=asin(sin(theta1)/surf->ior);

				VectAddS(-Refract_factor,inter->N,Zero_vect,L_vert);
				VectAddS(1/cos(theta2),L_vert,Zero_vect,Refract_vect);

				VectCopy(inter_ray->dir,Refract_vect);
				VectCopy(inter_ray->orig,Modify_P);

				trace_ray(level+1,weight,inter_ray,color);

		}
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// ray trace another pixel if the image isn't finished yet

void idle()
{
	if (image_j < ray_cam->im->h) {

		raytrace_one_pixel(image_i, image_j);

		image_i++;

		if (image_i == ray_cam->im->w) {
			image_i = 0;
			image_j++;
		}    
	}

	// write rendered image to file when done

	else if (!wrote_image) {

		write_PPM("output.ppm", ray_cam->im);

		wrote_image = true;
	}

	glutPostRedisplay();
}

//----------------------------------------------------------------------------

// show the image so far

void display(void)
{
	// draw it!

	glPixelZoom(1, -1);
	glRasterPos2i(0, ray_cam->im->h);

	glDrawPixels(ray_cam->im->w, ray_cam->im->h, GL_RGBA, GL_FLOAT, ray_cam->im->data);

	glFlush ();
}

//----------------------------------------------------------------------------

void init()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, ray_cam->im->w, 0.0, ray_cam->im->h);
}

//----------------------------------------------------------------------------

int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	// initialize scene (must be done before scene file is parsed)

	init_raytracing();

	if (argc == 2)
		parse_scene_file(argv[1], ray_cam);
	else {
		printf("missing .scene file\n");
		exit(1);
	}

	// opengl business

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ray_cam->im->w, ray_cam->im->h);
	glutInitWindowPosition(500, 300);
	glutCreateWindow("raytracing");
	init();

	glutDisplayFunc(display); 
	glutIdleFunc(idle);

	glutMainLoop();

	return 0; 
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
