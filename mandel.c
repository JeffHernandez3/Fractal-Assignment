
#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int threads);
void * run_me(void *arg);

struct parameters{
	int thread_id;
	int height;
	int width;
	int wide;
	double xmax;
	double xmin;
	double ymax;
	double ymin;
	int max;
	int threadMany;	//How many threads there are
	struct bitmap *part;
}params;

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main( int argc, char *argv[] )
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;

	//Number of threads used, where 1 is the default
	int threadNum = 1;
	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:n:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'n':
				threadNum = atoi(optarg);		//Changes the number of threads
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s threads= %d\n",xcenter,ycenter,scale,max,outfile,threadNum);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	
	// Compute the Mandelbrot image, this is the default if -n is not called
	if(threadNum==1)
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,threadNum);
	//start of threading 
	else if(threadNum > 0){
		//Makes spacing but more creating the parameters
		pthread_t tid[threadNum];
		struct parameters params[threadNum];

		int i;
		//Input information
		for(int i = 0; i < threadNum; i++){
			//struct bitmap *part = bitmap_create(image_width,image_height);			//Remember to delete
			params[i].xmin = xcenter - scale;
			params[i].xmax = xcenter + scale;
			params[i].ymin = ycenter - scale;
			params[i].ymax = ycenter + scale;
			params[i].max = max;
			params[i].threadMany = threadNum;
			params[i].thread_id = i;
			params[i].height = image_height;
			params[i].width = image_width;
			params[i].part = bm;
			pthread_create( &tid[i], NULL, run_me, (void*) &params[i]);
		}
		for (i = 0; i <threadNum; i++){
			pthread_join(tid[i], NULL);
		}
	}

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int threads)
{
	int i,j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	// For every pixel in the image...

	for(j=0;j<height;j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;
			//printf("\nThe xmin is %f, the j is %d, and the width is %d and finally the xmax is %f so my x is %f", xmin, j, width,xmax,x);

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}
	
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}
//Allows the thread to run 
void * run_me(void *arg){
	//Calling all of the parameters that seems to be needed
	int i;	//Checks for height
	int j;	//Checks for width
	double xmin;
	double xmax;
	double ymax;
	double ymin;
	int max;
	int threads;
	int height;
	int width;
	struct parameters *params = (struct parameters *) arg;
	//Putting values to those that were called beforehand
	xmin = params -> xmin;
	xmax = params -> xmax;
	ymax = params ->ymax;
	ymin = params ->ymin;
	max = params->max;
	threads = params ->threadMany;
	height = params ->height;
	width = params ->width;

	//Sets up how many threads we are going to use (height and width)
	//We need to get each pixel spot to be able to be used to then use iterations_at_point (such as 0 - 99 and such)
	int begin = params->thread_id *height / threads;
	int end = (begin + height/threads) + 1;

	//Starts the making of the image (similar to compute_image)
	//Difference is that j and i are flipped
	for(i = begin; i < end; i++){
		for(j = 0; j < width; j++){
			//Gets x and y at that spot
			double x = xmin + j * (xmax-xmin)/width;
			double y = ymin + i * (ymax-ymin)/height;
			//printf("\nThe xmin is %d, the j is %d, and the width is %d and finally the xmax is %d so my x is %f", xmin, j, width,xmax,x);
			//Compute iterations at that point
			int iters = iterations_at_point(x,y,max);
			//printf("\nChecking points for x: %f and y: %f", x, y);
			// Set the pixel in the bitmap of that thread
			bitmap_set(params->part,j,i,iters);
		}
		
	}

	return NULL;

}

