/*
Author: Liyang@Tsinghua
Date: 2018.12.2
Description: Automatically product KPOINTS file use POSCAR
Usage: gcc -o kpgen.c kpgen -lm;
       ./kpgen -s [kpoints_separation] (-f); 
*/

#include <stdio.h>
#include <stdlib.h>    
#include <string.h>    // Support: strchr
#include <unistd.h>    // Support: optget
#include <math.h>      // Support: sqrt, abs
                       // Add -lm if use 'gcc' compile this code

//##########################//
//### Function Statement ###//
//##########################//

// Function calculate the crossing multiply
void _CalCrossMultiply(double *, double *, double *);
// Function calculate the volume of a cell
double CalCellVolume(double *, double *, double *);
// Function calculate the length of a single reciprocal vector
double CalRecipVectorLength(double, double *);
// Function calculate the kpoints quantities in a reciprocal vector direction
int CalKpointsQuantity(double, double);

//#####################//
//### Main Function ###//
//#####################//

int main(int argc , char * argv[]){
  // Basic variable statement 
  double real_vector_a[3];
  double real_vector_b[3];
  double real_vector_c[3];
  double scaling_factor;
  double cross_real_vector_bc[3];
  double cross_real_vector_ca[3];
  double cross_real_vector_ab[3];
  double recip_vector_a[3];
  double recip_vector_b[3];
  double recip_vector_c[3];
  double length_recip_vector_a;
  double length_recip_vector_b;
  double length_recip_vector_c;
  double cell_volume;
  double k_points_separation = 0.03;
  int    kpoints_quantity_in_ka;
  int    kpoints_quantity_in_kb;
  int    kpoints_quantity_in_kc;
  int    is_2d_material = 0;

  // Read in the K points separation
  char *optstr = "s:f";
  int  opt;
  while((opt = getopt(argc, argv, optstr)) != -1){
    switch(opt){
      case 's':
        k_points_separation = atof(optarg);
        break;
      case 'f':
        is_2d_material = 1; 
        break;
      case '?':
        if(strchr(optstr, optopt) == NULL){
          fprintf(stderr, "[error] Unknown option '-%c'!!!\n", optopt);
        }else{
          fprintf(stderr, "[error] Option '-%c' requires an argument!!! \n", \
                  optopt);
        }
        return 0;
        }
    }

  // Read the info of real vector from the POSCAR 
  FILE *fpr = NULL;
  char buff[1024];
  if (fpr = fopen("POSCAR", "r")){
    fgets(buff, 1024, (FILE*)fpr);  // Skip the first line
    fgets(buff, 1024, (FILE*)fpr);  // Read the scaling factor
    scaling_factor = atof(buff);
    fscanf(fpr, "%s", buff);        // Read a_x
    real_vector_a[0] = atof(buff) * scaling_factor;
    fscanf(fpr, "%s", buff);        // Read a_y
    real_vector_a[1] = atof(buff) * scaling_factor;
    fscanf(fpr, "%s", buff);        // Read a_z
    real_vector_a[2] = atof(buff) * scaling_factor; 
    fscanf(fpr, "%s", buff);        // Read b_x
    real_vector_b[0] = atof(buff) * scaling_factor;
    fscanf(fpr, "%s", buff);        // Read b_y
    real_vector_b[1] = atof(buff) * scaling_factor;
    fscanf(fpr, "%s", buff);        // Read b_z
    real_vector_b[2] = atof(buff) * scaling_factor;
    fscanf(fpr, "%s", buff);        // Read c_x
    real_vector_c[0] = atof(buff) * scaling_factor;
    fscanf(fpr, "%s", buff);        // Read c_y
    real_vector_c[1] = atof(buff) * scaling_factor;
    fscanf(fpr, "%s", buff);        // Read c_z
    real_vector_c[2] = atof(buff) * scaling_factor;
    fclose(fpr);
  }else{
    fprintf(stderr, "[error] POSCAR do not exist in current folder!!!\n");
    return 1;
  }
  // Calculate the cell volume 
  cell_volume = CalCellVolume(real_vector_a, real_vector_b, real_vector_c);

  // Calculate the length of a reciprocal vector
  // -- Calculate the length of k_a
  _CalCrossMultiply(real_vector_b, real_vector_c, cross_real_vector_bc);
  length_recip_vector_a = CalRecipVectorLength(cell_volume, \
                                               cross_real_vector_bc);
  // -- Calculate the length of k_b
  _CalCrossMultiply(real_vector_c, real_vector_a, cross_real_vector_ca);
  length_recip_vector_b = CalRecipVectorLength(cell_volume, \
                                               cross_real_vector_ca);
  // -- Calculate the length of k_c
  _CalCrossMultiply(real_vector_a, real_vector_b, cross_real_vector_ab);
  length_recip_vector_c = CalRecipVectorLength(cell_volume, \
                                               cross_real_vector_ab);
                                               
  // Caculate the kpoints quantity in each reciprocal vector direction
  kpoints_quantity_in_ka = CalKpointsQuantity(length_recip_vector_a,\
                                              k_points_separation);
  kpoints_quantity_in_kb = CalKpointsQuantity(length_recip_vector_b,\
                                              k_points_separation);
  if(is_2d_material == 0){ 
    kpoints_quantity_in_kc = CalKpointsQuantity(length_recip_vector_c,\
                                                k_points_separation);
  }else{
    kpoints_quantity_in_kc = 1;
  }

  // Write the KPOINTS file
  FILE * fpw = NULL;
  fpw = fopen("KPOINTS", "w");
  fprintf(fpw, "Auto.Mesh.Kpts\n");
  fprintf(fpw, "0             \n");
  fprintf(fpw, "Monkhorst-Pack\n");
  fprintf(fpw, "%d    %d    %d\n", kpoints_quantity_in_ka,\
                                   kpoints_quantity_in_kb,\
                                   kpoints_quantity_in_kc);
  fprintf(fpw, "0.0  0.0  0.0 \n");
  fclose(fpw);

  return 0;
}

//###########################//
//### Function Defination ###//
//###########################//

// Function calculate the crossing multiply
void _CalCrossMultiply(double *A, double *B, double *cross_multiply_A_B){
  cross_multiply_A_B[0] = A[1] * B[2] - A[2] * B[1];
  cross_multiply_A_B[1] = A[2] * B[0] - A[0] * B[2];
  cross_multiply_A_B[2] = A[0] * B[1] - A[1] * B[0];
}

// Function calculate the volume of a cell
double CalCellVolume(double *real_a, double *real_b, double *real_c){
  double cell_volume;
  double cross_multiply_bc[3];
  _CalCrossMultiply(real_b, real_c, cross_multiply_bc);
  cell_volume = real_a[0] * cross_multiply_bc[0] + \
                real_a[1] * cross_multiply_bc[1] + \
                real_a[2] * cross_multiply_bc[2];
  return cell_volume;
}

// Function calculate the length of a single reciprocal vector
double CalRecipVectorLength(double cell_volume, double *cross_multiply_jk){
  double recip_vector_i[3];
  double length_recip_vector_i;
  // Calculate each component of the reciprocal vector
  recip_vector_i[0] = (1.0 / cell_volume) * cross_multiply_jk[0];
  recip_vector_i[1] = (1.0 / cell_volume) * cross_multiply_jk[1];
  recip_vector_i[2] = (1.0 / cell_volume) * cross_multiply_jk[2];
  // Calculate the length of this reciprocal vector
  length_recip_vector_i = sqrt(recip_vector_i[0] * recip_vector_i[0] + \
                               recip_vector_i[1] * recip_vector_i[1] + \
                               recip_vector_i[2] * recip_vector_i[2]);
  // Return the absolute value of the calcualtion
  return fabs(length_recip_vector_i);
}

// Function calculate the kpoints quantities in a reciprocal vector direction
int CalKpointsQuantity(double length_recip_vector_i, \
                       double k_points_separation){
  int kpoints_quantity_in_ki;
  kpoints_quantity_in_ki = (int)(length_recip_vector_i / k_points_separation);
  if(kpoints_quantity_in_ki % 2 == 0){
    if(kpoints_quantity_in_ki == 2){
      kpoints_quantity_in_ki--;
    }else{
      kpoints_quantity_in_ki++;
    }
  }
  return kpoints_quantity_in_ki;
}

