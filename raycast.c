/*
 *raycast.c
 *By: Logan Green
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>


#define MAX_DEPTH 7   
//Creating the various parts of the scene
typedef struct{
  int type; 
  double position[3];
  double diffuse_color[3];
  double specular_color[3];
  double reflectivity	;
  double refractivity;
  double ior;
  union {
    struct {
      double radius;
    } sphere;
    struct {
      double normal[3];
    } plane;
  };
} Object;

typedef struct{
  double width;
  double height;
} Camera;

typedef struct{
  double color[3];
  double position[3];
  double direction[3];
  double radial_a0;
  double radial_a1;
  double radial_a2;
  double theta;
  double angular_a0;
} Light;

typedef struct Pixel{
  unsigned char r, b, g;
}Pixel;

typedef struct Closest{
	Object* closest_object;
	double closest_t;
}Closest;

//JSON reading functions
int next_c(FILE* json);
void check_c(FILE* json, int d);
void skip_ws(FILE* json);
char* next_string(FILE* json);
double next_number(FILE* json);
double* next_vector(FILE* json);
void read_scene(char* filename, Camera* camera, Object** objects, Light** lights);

//Functions for vectors
void vector_normalize(double* v);
double vector_dot_product(double *v1, double *v2);
void vector_cross_product(double *v1, double *v2, double *result);
double vector_length(double *vector);
void vector_reflection(double *N, double *L, double *result);
void vector_subtraction(double *v1, double *v2, double *result);
void vector_addition(double *v1, double *v2, double *result);
void vector_scale(double *vector, double scalar, double *result);

//Functions for intersections
double sphere_intersection(double* Ro, double* Rd, double* C, double r);
double plane_intersection(double* Ro, double* Rd, double* P, double* N);
Closest* shoot(double* Ro, double* Rd, Object **objects);

//Functions for light
double calculate_diffuse(double object_diff_color, double light_color, double *N, double *L);
double calculate_specular(double *L, double *N, double *R, double *V, double object_spec_color, double light_color);
double frad(Light * light, double t);
double fang(Light *light, double *L);

//Functions for images
void generate_scene(Camera* camera, Object** objects, Light** lights, Pixel* buffer, int width, int height);
Pixel* recursive_shade(Object **objects, Light **lights, double* Ro, double* Rd, Closest* current_object, int depth, double current_ior, int exiting_sphere);
void write_p3(Pixel *buffer, FILE *output_file, int width, int height, int max_color);
double clamp(double value);

//Where they at tho.
int line = 1;

int main(int argc, char *argv[]) {
  //Checks to see if there are five arguments
  if (argc != 5){
    fprintf(stderr, "Error: 5 arguments expected. there are currently: %d.\n", argc);
  }
  int width = atoi(argv[1]);
  int height = atoi(argv[2]);
  if (width <= 0 || height <=){
    fprintf(stderr, "Error: height or width is not in the proper dimentions.\n");
  }
  Object **objects;
  objects = malloc(sizeof(Object*)*128);
  //create camera object
  Camera *camera;
  camera = (Camera *)malloc(sizeof(Camera));
  //create lights
  Light **lights;
  lights = malloc(sizeof(Light*)*128);
  //create image buffer
  Pixel *buffer; 
  buffer = (Pixel *)malloc(width*height*sizeof(Pixel));
  read_scene(argv[3], camera, objects, lights);
  generate_scene(camera, objects, lights, buffer, width, height);
  FILE* output_file = fopen(argv[4], "w");
  if (output_file == NULL){
    fprintf(stderr, "Error: Unable to open file.\n");
    fclose(output_file);
    exit(1);
  }
  write_p3(buffer, output_file, width, height, 255);
  fclose(output_file);
  for(int i = 0; lights[i]!=NULL; i+=1){
  	free(objects[i]);
  }
  free(objects);
  for(int i = 0; lights[i]!=NULL; i+=1){
  	free(lights[i]);
  }
  free(lights);
  free(camera);
  return EXIT_SUCCESS;
}
//essentially the getc() function but with more usability
int next_c(FILE *json) {
  int c = fgetc(json);
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number: %d.\n", line);
    fclose(json);
    exit(1);
  }
  return c;
}

//checks that the next character is correct
void check_c(FILE *json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error.\n");
  fclose(json);
  exit(1);    
}

//simple function to skip whitespace
void skip_ws(FILE *json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}
//gets the next string
char* next_string(FILE *json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line: %d.\n", line);
    fclose(json);
    exit(1);
  }  
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: String length is too long.\n");
      fclose(json);
      exit(1);      
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE *json) {
  double value;
  int count = fscanf(json, "%lf", &value);
  if (count != 1){
    fprintf(stderr, "Error: Failed to read number.\n");
    fclose(json);
    exit(1);
  }
  return value;
}

double* next_vector(FILE *json) {
  double* v = malloc(3*sizeof(double));
  check_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  check_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  check_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  check_c(json, ']');
  return v;
}

void read_scene(char *filename, Camera *camera, Object **objects, Light **lights) {
  int c;
  int current_light = -1;
  int current_item = -1;
  int current_type;
  FILE* json = fopen(filename, "r");
  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file.\n");
    fclose(json);
    exit(1);
  }

  skip_ws(json);
  check_c(json, '[');
  skip_ws(json);

  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: The scene file is too short.\n");
      fclose(json);
      return;
    }
    if (c == '{') {
      skip_ws(json);
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
        fprintf(stderr, "Error: Expected type on line number %d.\n", line);
        exit(1);
      }

      skip_ws(json);
      check_c(json, ':');
      skip_ws(json);

      char* value = next_string(json);
      if(strcmp(value, "camera") == 0){
        current_type = 0;
      } 
      else if(strcmp(value, "sphere") == 0) {
        current_item++;
        if(current_item>127){
          fprintf(stderr, "Error: Too many objects.\n");
          fclose(json);
          exit(1);
        }
        objects[current_item] = malloc(sizeof(Object));
        objects[current_item]->type = 0;
        objects[current_item]->reflectivity = 0.0;
        objects[current_item]->refractivity = 0.0;
        objects[current_item]->ior = 1.0;
        current_type = 1;
      } 
      else if (strcmp(value, "plane") == 0) {
        current_item++;
        if(current_item>127){
          fprintf(stderr, "Error: Too many objects.\n");
          fclose(json);
          exit(1);
        }
        objects[current_item] = malloc(sizeof(Object));
        objects[current_item]->type = 1;
        objects[current_item]->reflectivity = 0.0;
        objects[current_item]->refractivity = 0.0;
        objects[current_item]->ior = 1.0;
        current_type = 2;
      }
      else if (strcmp(value, "light") == 0) {
        current_light++;
        if(current_light>127){
          fprintf(stderr, "Error: Too many lights.\n");
          fclose(json);
          exit(1);
        }
        lights[current_light] = malloc(sizeof(Light));
        current_type = 3;
      } 
      else { 
        fprintf(stderr, "Error: Unknown type.\n");
        fclose(json);
        exit(1);
      }

      skip_ws(json);

      while (1) {
        c = next_c(json);
        if (c == '}') {
          break;
        } 
        else if (c == ',') {
          skip_ws(json);
          char* key = next_string(json);
          skip_ws(json);
          check_c(json, ':');
          skip_ws(json);
          if (strcmp(key, "width") == 0){
            if(current_type == 0){
              camera->width = next_number(json);
            }else{
              fprintf(stderr, "Error: Current width value is not correct.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "height") == 0){
            if(current_type == 0){
              camera->  height = next_number(json);
            }else{
              fprintf(stderr, "Error: Current height value is not correct.\n");
              fclose(json);
              exit(1);
            }
          }else if (strcmp(key, "radial-a0") == 0){
            if(current_type == 3){
              lights[current_light]->radial_a0 = next_number(json);
            }else{
              fprintf(stderr, "Error\n");
              fclose(json);
              exit(1);
            }
          }else if (strcmp(key, "radial-a1") == 0){
            if(current_type == 3){
              lights[current_light]->radial_a1 = next_number(json);
            }else{
              fprintf(stderr, "Error\n");
              fclose(json);
              exit(1);
            }
          }else if (strcmp(key, "radial-a2") == 0){
            if(current_type == 3){
              lights[current_light]->radial_a2 = next_number(json);
            }else{
              fprintf(stderr, "Error\n");
              fclose(json);
              exit(1);
            }
          }else if (strcmp(key, "angular-a0") == 0){
            if(current_type == 3){
              lights[current_light]->angular_a0 = next_number(json);
            }else{
              fprintf(stderr, "Error\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "radius") == 0){
            if(current_type == 1){
              objects[current_item]->sphere.radius = next_number(json);
            }else{
              fprintf(stderr, "Error: Current object cannot have radius value.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "diffuse_color") == 0){ 
            if(current_type == 1 || current_type == 2){
                double* vector = next_vector(json);
                objects[current_item]->diffuse_color[0] = vector[0];
                objects[current_item]->diffuse_color[1] = vector[1];
                objects[current_item]->diffuse_color[2] = vector[2];
                free(vector);
            }else{
              fprintf(stderr, "Error: Non-object type has color value.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "specular_color") == 0){ 
            if(current_type == 1 || current_type == 2){
                double* vector = next_vector(json);
                objects[current_item]->specular_color[0] = vector[0];
                objects[current_item]->specular_color[1] = vector[1];
                objects[current_item]->specular_color[2] = vector[2];
                free(vector);
            }else{
              fprintf(stderr, "Error: Non-object type has color value.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "reflectivity") == 0){
            if(current_type == 1 || current_type == 2){
                objects[current_item]->reflectivity = next_number(json);
                int value = objects[current_item]->reflectivity + objects[current_item]->refractivity;
                if (value > 1){
                	fprintf(stderr, "Error: reflectivity is too high\n");
              		fclose(json);
              		exit(1);
                }
            }else{
              fprintf(stderr, "Error: Non-object type has reflectivity value.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "refractivity") == 0){ 
            if(current_type == 1 || current_type == 2){
                objects[current_item]->refractivity = next_number(json);
                int value = objects[current_item]->reflectivity + objects[current_item]->refractivity;
                if (value > 1){
                	fprintf(stderr, "Error: Refractivity is to high.\n");
              		fclose(json);
              		exit(1);
                }
            }else{
              fprintf(stderr, "Error: Non-object type has refractivity value on line number %d.\n", line);
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "ior") == 0){ 
            if(current_type == 1 || current_type == 2){
                objects[current_item]->ior = next_number(json);
            }else{
              fprintf(stderr, "Error: Non-object type has IoR value.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "color") == 0){ 
            if(current_type == 3){
                double* vector = next_vector(json);
                lights[current_light]->color[0] = vector[0];
                lights[current_light]->color[1] = vector[1];
                lights[current_light]->color[2] = vector[2]; 
                free(vector);
            }else{
              fprintf(stderr, "Error: Non-light type has color value.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "position") == 0){
            if(current_type == 1 || current_type == 2){
              double* vector = next_vector(json);
              objects[current_item]->position[0] = vector[0];
              objects[current_item]->position[1] = vector[1];
              objects[current_item]->position[2] = vector[2];
              free(vector);
            }else if(current_type == 3){
              double* vector = next_vector(json);
              lights[current_light]->position[0] = vector[0];
              lights[current_light]->position[1] = vector[1];
              lights[current_light]->position[2] = vector[2];  
              free(vector);
            }else{
              fprintf(stderr, "Error.\n");
              fclose(json);
              exit(1);
            }
          } else if(strcmp(key, "normal") == 0){
            if(current_type == 2){
              double* vector = next_vector(json);
              objects[current_item]->plane.normal[0] = vector[0];
              objects[current_item]->plane.normal[1] = vector[1];
              objects[current_item]->plane.normal[2] = vector[2];  
              free(vector);
            }else{
              fprintf(stderr, "Error: Non-Plane has normal value.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "direction") == 0){
            if(current_type == 3){
              double* vector = next_vector(json);
              lights[current_light]->direction[0] = vector[0];
              lights[current_light]->direction[1] = vector[1];
              lights[current_light]->direction[2] = vector[2];  
              free(vector);
            }else{
              fprintf(stderr, "Error: Non-Light has direction.\n");
              fclose(json);
              exit(1);
            }
          }else if(strcmp(key, "theta") == 0){
            if(current_type == 3){
              lights[current_light]->theta = next_number(json);
            }else{
              fprintf(stderr, "Error: Non-Sphere has a theta value.\n");
              fclose(json);
              exit(1);
            }
          }else{
            fprintf(stderr, "Error: Unknown value.\n",
                key, line);
            fclose(json);
            exit(1);
          }

          skip_ws(json);
        }else {
          fprintf(stderr, "Error: Unexpected value.\n");
          fclose(json);
          exit(1);
        }
      }

      skip_ws(json);
      c = next_c(json);
      if (c == ',') {
        skip_ws(json);
      } 
      else if (c == ']') {
        fclose(json);
        return;
      } 
      else {
        fprintf(stderr, "Error: Expected ',' or ']'.\n");
        fclose(json);
        exit(1);
      }
    }
  }
}

void vector_normalize(double *v) {
  double len = sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

double vector_dot_product(double *v1, double *v2){
	return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

void vector_cross_product(double *v1, double *v2, double *result){
	result[0] = v1[1]*v2[2]-v2[1]*v1[2];
	result[1] = v2[0]*v1[2]-v1[0]*v2[2];
	result[2] = v1[0]*v2[1]-v2[0]*v1[1];
}

double vector_length(double *vector){
	return sqrt(pow(vector[0], 2)+pow(vector[1], 2)+pow(vector[2], 2));
}

void vector_reflection(double *N, double *L, double *result){
	double dot_result;
	double temp_vector[3]; 
	dot_result = 2*vector_dot_product(N, L);
	vector_scale(N, dot_result, temp_vector);
	vector_subtraction(L, temp_vector, result);

}

void vector_subtraction(double *v1, double *v2, double *result){
	result[0] = v2[0] - v1[0];
	result[1] = v2[1] - v1[1];
	result[2] = v2[2] - v1[2];
}

void vector_addition(double *v1, double *v2, double *result){
	result[0] = v2[0] + v1[0];
	result[1] = v2[1] + v1[1];
	result[2] = v2[2] + v1[2];
}


void vector_scale(double *vector, double scalar, double *result){
	result[0] = vector[0]*scalar;
	result[1] = vector[1]*scalar;
	result[2] = vector[2]*scalar;
}

double sphere_intersection(double *Ro, double *Rd, double *C, double r) {
  double a = (pow(Rd[0], 2) + pow(Rd[1], 2) + pow(Rd[2], 2));
  double b = (2 * (Ro[0] * Rd[0] - Rd[0] * C[0] + Ro[1] * Rd[1] - Rd[1] * C[1] + Ro[2] * Rd[2] - Rd[2] * C[2]));
  double c = pow(Ro[0], 2) - 2*Ro[0]*C[0] + pow(C[0], 2) + pow(Ro[1], 2) - 2*Ro[1]*C[1] + pow(C[1], 2) + pow(Ro[2], 2) - 2*Ro[2]*C[2] + pow(C[2], 2) - pow(r, 2);

  double det = pow(b, 2) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);

  double t0 = (-b - det) / (2*a);
  if (t0 > 0.00001) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0.00001) return t1;

  return -1;
}

double plane_intersection(double *Ro, double *Rd, double *P, double *N) {
  vector_normalize(N);
  double t = (N[0]*P[0] + N[1]*P[1] + N[2]*P[2] - N[0]*Ro[0] - N[1]*Ro[1] - N[2]*Ro[2])/(N[0]*Rd[0] + N[1]*Rd[1] + N[2]*Rd[2]); 
  if (t > 0) return t;

  return -1;
}

Closest* shoot(double *Ro, double *Rd, Object **objects){
	Closest* best_values = malloc(sizeof(Closest));
	best_values->closest_object = NULL;
	best_values->closest_t = INFINITY;
	vector_normalize(Rd);
  	for (int i=0; objects[i] != 0; i += 1) {
		double t = 0;
		switch(objects[i]->type) {
			case 0:
		  		t = sphere_intersection(Ro, Rd, objects[i]->position, objects[i]->sphere.radius);
		  		break;
			case 1:
		  		t = plane_intersection(Ro, Rd, objects[i]->position, objects[i]->plane.normal);
		  		break;
			default:
		  		printf("Error: Unknown object type.\n");	
		  		exit(1);
		}
		if (t > 0.00001 && t < best_values->closest_t) {
			best_values->closest_t = t;
	  		best_values->closest_object = objects[i];
		}
	}
	return best_values;
}

double calculate_diffuse(double object_diff_color, double light_color, double *N, double *L){
	double dot_result; 
	dot_result = vector_dot_product(N, L);
	if(dot_result > 0){
		return object_diff_color*light_color*dot_result;
	}
	else{
		return 0;
	}
}

double calculate_specular(double *L, double *N, double *R, double *V, double object_spec_color, double light_color){
	double V_dot_R;
	double N_dot_L;
	V_dot_R = vector_dot_product(V, R);
	N_dot_L = vector_dot_product(N, L);
	if(V_dot_R > 0 && N_dot_L > 0){
		return object_spec_color * light_color * pow(V_dot_R, 20);
	}
	else{
		return 0;
	}
}

double frad(Light * light, double t){
	return (1/(light->radial_a2*t*t+light->radial_a1*t+light->radial_a0));
}

double fang(Light *light, double *L){
	double light_length = vector_length(light->direction);
	double theta = light->theta;
	double dot_result;
	double light_vector[3];

	theta = theta*M_PI/180;
	theta = cos(theta); 
	
	light_vector[0] = light->direction[0];
	light_vector[1] = light->direction[1];
	light_vector[2] = light->direction[2];
	
	vector_normalize(light_vector);
	
	if(light->theta == 0 || light_length == 0){
		return 1;
	}
	else{	
		dot_result = vector_dot_product(light_vector, L);
		if (dot_result > theta){
			return 0;
		}
		else{
			return pow(dot_result, light->angular_a0);
		}
	}
	return 0;
}

void generate_scene(Camera *camera, Object **objects, Light **lights, Pixel *buffer, int width, int height){
  double camera_width = camera->width;
  double camera_height = camera->height;
  double pixheight = camera_height / height;
  double pixwidth = camera_width / width;
  Pixel* current_pixel;
  int position;
  for (int y = 0; y < height; y += 1) {
    for (int x = 0; x < width; x += 1) {
      	double Ro[3] = {0, 0, 0};
      	double Rd[3] = {
        	0 - (camera_width/2) + pixwidth * (x + 0.5),
        	0 - (camera_height/2) + pixheight * (y + 0.5),
        	1
  		};
      	vector_normalize(Rd);
  		Closest* nearest_object = shoot(Ro, Rd, objects);
		if (nearest_object->closest_t > 0 && nearest_object->closest_t != INFINITY) {
			current_pixel = recursive_shade(objects, lights, Ro, Rd, nearest_object, 0, 1.0, 0);
		}	 
		else {
		  	current_pixel->r = 0;
		  	current_pixel->g = 0;
		  	current_pixel->b = 0;
		}
      	position = (height-(y+1))*width+x;
      	buffer[position].r = current_pixel->r;
      	buffer[position].g = current_pixel->g;
      	buffer[position].b = current_pixel->b;
      	free(current_pixel);
    }
  } 
}

Pixel* recursive_shade(Object **objects, Light **lights, double *Ro, double *Rd, Closest *current_object, int depth, double current_ior, int exiting_sphere){
	Pixel* current_pixel = malloc(sizeof(Pixel));
	Object* closest_object = current_object->closest_object;
	double closest_t = current_object->closest_t;
	double color[3];
	double Ron[3];
	double Rdn[3];
	double N[3];
	double L[3];
	double R[3];
	double V[3];
	double radial_light;
	double angular_light;
	color[0] = 0;
  	color[1] = 0;
  	color[2] = 0;
  	Pixel* reflect = malloc(sizeof(Pixel));
  	reflect->r = 0;
  	reflect->g = 0;
  	reflect->b = 0;
  	Pixel* refract = malloc(sizeof(Pixel));
  	refract->r = 0;
  	refract->g = 0;
  	refract->b = 0;
  		
  	if(current_object->closest_object->reflectivity > 0.00001 && depth <= MAX_DEPTH){
  	  double new_ray[3];
  	  vector_scale(Rd, -1, new_ray);
  	  Ron[0] = closest_t * Rd[0] + Ro[0];
      	  Ron[1] = closest_t * Rd[1] + Ro[1];
      	  Ron[2] = closest_t * Rd[2] + Ro[2];
      	  if(closest_object->type == 1){
	    N[0] = closest_object->plane.normal[0];
	    N[1] = closest_object->plane.normal[1];
	    N[2] = closest_object->plane.normal[2];
         }else if(closest_object->type == 0){
	    N[0] = Ron[0] - closest_object->position[0];
	    N[1] = Ron[1] - closest_object->position[1];
	    N[2] = Ron[2] - closest_object->position[2];
         }else{
	    printf("Error: Unknown object type.\n");	
            exit(1);
         }
	vector_normalize(N);
     	vector_normalize(new_ray);
       	vector_reflection(N, new_ray, R);
       	vector_normalize(R);
	Closest* next_surface = shoot(Ron, R, objects);	
	if(next_surface->closest_t > 0 && next_surface->closest_t < INFINITY){;
          int new_depth = depth+1;
      	  reflect = recursive_shade(objects, lights, Ron, R, next_surface, new_depth, current_ior, 0);
       	}
      }
      if(closest_object->refractivity > 0.00001 && depth <= MAX_DEPTH){
  	double new_origin[3];
  	new_origin[0] = closest_t * Rd[0] + Ro[0];
      	new_origin[1] = closest_t * Rd[1] + Ro[1];
      	new_origin[2] = closest_t * Rd[2] + Ro[2];
       	double new_ray[3] = {Rd[0], Rd[1], Rd[2]};
	double a[3];
       	double b[3];
       	double sin_theta;
       	double sin_phi;
       	double cos_phi;
       	double external_ior;
       	double ior;

       	if(closest_object->type == 1){
	  N[0] = closest_object->plane.normal[0];
	  N[1] = closest_object->plane.normal[1];
	  N[2] = closest_object->plane.normal[2];
        }else if(closest_object->type == 0){
	  N[0] = Ro[0] - closest_object->position[0];
	  N[1] = Ro[1] - closest_object->position[1];
	  N[2] = Ro[2] - closest_object->position[2];
	}else{
	  printf("Error: Unknown object type.\n");	
          exit(1);
	}
	if(exiting_sphere == 1){
	  external_ior = closest_object->ior*current_ior;
	  ior = closest_object->ior/external_ior;
	}else{
	  ior = current_ior/closest_object->ior;
	}
       	vector_normalize(N);
       	vector_normalize(Rd);
       	vector_cross_product(N, Rd, a);
       	vector_normalize(a);
       	vector_cross_product(a, N, b);
       	vector_normalize(b);
       	sin_theta = vector_dot_product(Rd, b);
       	sin_phi = ior*sin_theta;
       	if(pow(sin_phi, 2)<=1){
          cos_phi = sqrt(1-pow(sin_phi, 2));	  		
	  vector_scale(N, -1*cos_phi, N);
	  vector_scale(b, sin_phi, b);
	  vector_addition(N, b, new_ray);
	  vector_normalize(new_ray);
	  Closest* next_surface = shoot(new_origin, new_ray, objects);
	  if(next_surface->closest_t > 0 && next_surface->closest_t < INFINITY){
	    int new_depth = depth + 1;
	    if(next_surface->closest_object == closest_object){
	      refract = recursive_shade(objects, lights, new_origin, new_ray, next_surface, new_depth, ior, 1);
	    }else{
	      if(exiting_sphere == 1){
		refract = recursive_shade(objects, lights, new_origin, new_ray, next_surface, new_depth, external_ior, 0);
	      }else{
		refract = recursive_shade(objects, lights, new_origin, new_ray, next_surface, new_depth, ior, 0);
	      }
	    }
	  }
	}
      }
      for (int j=0; lights[j] != NULL; j+=1){
        Ron[0] = closest_t * Rd[0] + Ro[0];
      	Ron[1] = closest_t * Rd[1] + Ro[1];
      	Ron[2] = closest_t * Rd[2] + Ro[2];
      	Rdn[0] = lights[j]->position[0] - Ron[0];
      	Rdn[1] = lights[j]->position[1] - Ron[1];
      	Rdn[2] = lights[j]->position[2] - Ron[2];
      	double distance_to_light = vector_length(Rdn);
      	vector_normalize(Rdn);
      	double shadow_t;
      	int closest_shadow_object = 0;
	if(closest_object->type == 1){
	  N[0] = closest_object->plane.normal[0];
	  N[1] = closest_object->plane.normal[1];
	  N[2] = closest_object->plane.normal[2];
 	}else if(closest_object->type == 0){
	  N[0] = Ron[0] - closest_object->position[0];
	  N[1] = Ron[1] - closest_object->position[1];
	  N[2] = Ron[2] - closest_object->position[2];
       	}else{
	  printf("Error: Unknown object type.\n");
	  exit(1);
       	}
	vector_normalize(N);

      	for (int k=0; objects[k] != NULL; k+=1) {
	  if (objects[k] == closest_object) {
	    continue;
	  }
	  switch(objects[k]->type) {
	  case 0:
	    shadow_t = sphere_intersection(Ron, Rdn, objects[k]->position, objects[k]->sphere.radius);
	    break;
	  case 1:
	    shadow_t = plane_intersection(Ron, Rdn, objects[k]->position, objects[k]->plane.normal);
	    break;
	  default:
	    printf("Error: Unknown object type.\n");
	    exit(1);
	  }
	  if (0 < shadow_t && shadow_t < distance_to_light) {
	    closest_shadow_object = 1;
	    break;
	  }
	}
	if (closest_shadow_object == 0) {
	  if(closest_object->type == 1){
	    N[0] = closest_object->plane.normal[0];
	    N[1] = closest_object->plane.normal[1];
	    N[2] = closest_object->plane.normal[2];
	  }else if(closest_object->type == 0){
	    N[0] = Ron[0] - closest_object->position[0];
	    N[1] = Ron[1] - closest_object->position[1];
	    N[2] = Ron[2] - closest_object->position[2];
	  }
	  else{
	    printf("Error: Unknown object type.\n");
	    exit(1);
	  }
	  vector_normalize(N);
	  L[0] = Rdn[0];
	  L[1] = Rdn[1];
	  L[2] = Rdn[2];
	  vector_normalize(L);
	  //Get R
	  vector_reflection(N, L, R);
	  vector_normalize(R);
	  //Get V
	  V[0] = -1*Rd[0];
	  V[1] = -1*Rd[1];
	  V[2] = -1*Rd[2];
	  vector_normalize(V);

	  double diffuse[3];
	  diffuse[0] = calculate_diffuse(closest_object->diffuse_color[0], lights[j]->color[0], N, L);
	  diffuse[1] = calculate_diffuse(closest_object->diffuse_color[1], lights[j]->color[1], N, L);
	  diffuse[2] = calculate_diffuse(closest_object->diffuse_color[2], lights[j]->color[2], N, L);

	  double specular[3];
	  specular[0] = calculate_specular(L, N, R, V, closest_object->specular_color[0], lights[j]->color[0]);
	  specular[1] = calculate_specular(L, N, R, V, closest_object->specular_color[1], lights[j]->color[1]);
	  specular[2] = calculate_specular(L, N, R, V, closest_object->specular_color[2], lights[j]->color[2]);
	  radial_light = frad(lights[j], distance_to_light);
	  angular_light = fang(lights[j], L);
	  color[0] += (radial_light * angular_light * (diffuse[0] + specular[0]));
	  color[1] += (radial_light * angular_light * (diffuse[1] + specular[1]));
	  color[2] += (radial_light * angular_light * (diffuse[2] + specular[2]));
      	}
      }
      double reflective[3];
      reflective[0] = ((double)reflect->r)/255;
      reflective[1] = ((double)reflect->g)/255;
      reflective[2] = ((double)reflect->b)/255;
      double refractive[3];
      refractive[0] = ((double)refract->r)/255;
      refractive[1] = ((double)refract->g)/255;
      refractive[2] = ((double)refract->b)/255;
      color[0] = (color[0])*(1-closest_object->reflectivity-closest_object->refractivity);
      color[0] += (closest_object->reflectivity*reflective[0]);
      color[0] += (closest_object->refractivity*refractive[0]);
      color[1] = (color[1])*(1-closest_object->reflectivity-closest_object->refractivity);
      color[1] += (closest_object->reflectivity*reflective[1]);
      color[1] += (closest_object->refractivity*refractive[1]);
      color[2] = (color[2])*(1-closest_object->reflectivity-closest_object->refractivity);
      color[2] += (closest_object->reflectivity*reflective[2]);
      color[2] += (closest_object->refractivity*refractive[2]);
      current_pixel->r = (unsigned char)(255 * clamp(color[0]));
      current_pixel->g = (unsigned char)(255 * clamp(color[1]));
      current_pixel->b = (unsigned char)(255 * clamp(color[2]));
      return current_pixel;
}

void write_p3(Pixel *buffer, FILE *output_file, int width, int height, int max_color){
  fprintf(output_file, "P3\n%d %d\n%d\n", width, height, max_color);
  int current_width = 1;
  for(int i = 0; i < width*height; i++){
    fprintf(output_file, "%d %d %d ", buffer[i].r, buffer[i].g, buffer[i].b);
    if(current_width >= 70%12){
      fprintf(output_file, "\n");
      current_width = 1;
    }
    else{
      current_width++;
    }
  }
}

double clamp(double value){
  if (value < 0){
    return 0;
  }else if (value > 1){
    return 1;
  }else{
    return value;
  }
}
