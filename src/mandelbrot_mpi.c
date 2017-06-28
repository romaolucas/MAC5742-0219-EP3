#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MASTER 0

double c_x_min;
double c_x_max;
double c_y_min;
double c_y_max;

double pixel_width;
double pixel_height;

int iteration_max = 200;

int image_size;
unsigned char **image_buffer;

int i_x_max;
int i_y_max;
int image_buffer_size;

int gradient_size = 16;
int colors[17][3] = {
                        {66, 30, 15},
                        {25, 7, 26},
                        {9, 1, 47},
                        {4, 4, 73},
                        {0, 7, 100},
                        {12, 44, 138},
                        {24, 82, 177},
                        {57, 125, 209},
                        {134, 181, 229},
                        {211, 236, 248},
                        {241, 233, 191},
                        {248, 201, 95},
                        {255, 170, 0},
                        {204, 128, 0},
                        {153, 87, 0},
                        {106, 52, 3},
                        {16, 16, 16},
                    };

void allocate_image_buffer(){
    int rgb_size = 3;
    image_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * image_buffer_size);

    for(int i = 0; i < image_buffer_size; i++){
        image_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
    };
};


void init(int argc, char *argv[]){
    if(argc < 6){
        printf("usage: mpirun -n num_processors mandelbrot_mpi c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500 and 10 processors:\n");
        printf("    Full Picture:        mpirun -n 10 mandelbrot_mpi -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:     mpirun -n 10 mandelbrot_mpi -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:     mpirun -n 10 mandelbrot_mpi 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: mpirun -n 10 mandelbrot_mpi -0.188 -0.012 0.554 0.754 11500\n");
        exit(0);
    }
    else{
        sscanf(argv[1], "%lf", &c_x_min);
        sscanf(argv[2], "%lf", &c_x_max);
        sscanf(argv[3], "%lf", &c_y_min);
        sscanf(argv[4], "%lf", &c_y_max);
        sscanf(argv[5], "%d", &image_size);
        i_x_max           = image_size;
        i_y_max           = image_size;
        image_buffer_size = image_size * image_size;

        pixel_width       = (c_x_max - c_x_min) / i_x_max;
        pixel_height      = (c_y_max - c_y_min) / i_y_max;
    };
};

void compute_mandelbrot(int y_min, int y_max){
    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;

    int iteration;
    int i_x;
    int i_y;

    double c_x;
    double c_y;
    for(i_y = y_min; i_y < y_max; i_y++){
        c_y = c_y_min + i_y * pixel_height;

        if(fabs(c_y) < pixel_height / 2){
            c_y = 0.0;
        };

        for(i_x = 0; i_x < i_x_max; i_x++){
            c_x         = c_x_min + i_x * pixel_width;

            z_x         = 0.0;
            z_y         = 0.0;

            z_x_squared = 0.0;
            z_y_squared = 0.0;

            for(iteration = 0;
                iteration < iteration_max && \
                ((z_x_squared + z_y_squared) < escape_radius_squared);
                iteration++){
                z_y         = 2 * z_x * z_y + c_y;
                z_x         = z_x_squared - z_y_squared + c_x;

                z_x_squared = z_x * z_x;
                z_y_squared = z_y * z_y;
            };

            //update_rgb_buffer(iteration, i_x, i_y);
        };
    };
};

int main(int argc, char *argv[]){
    int numtasks, taskid;
    int y_min, y_max;
    int step;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    
    init(argc, argv);
    if (taskid == MASTER) {
        //allocate_image_buffer();
        
        int i;
        step = image_size / numtasks; 
        for (i = 1; i < numtasks; i++) {
            //MPI_Send(&image_buffer, image_buffer_size*3, MPI_UNSIGNED_CHAR, i, MPI_COMM_WORLD);    
            y_min = i * step;
            y_max = y_min + step;
            if (i == numtasks - 1) 
                y_max += image_size % numtasks;
            MPI_Send(&y_min, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&y_max, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
        }
        y_min = 0;
        y_max = step;
        compute_mandelbrot(y_min, y_max);
    } else {
        //MPI_Recv(&image_buffer, image_buffer_size*3, MPI_UNSIGNED_CHAR, i, MPI_COMM_WORLD);
        MPI_Recv(&y_min, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&y_max, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        compute_mandelbrot(y_min, y_max);
    }
    //write_to_file();
    MPI_Finalize();

    return 0;
};
