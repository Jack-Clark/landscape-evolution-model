// This file contains compiler directives to allow the code to be compiled on differnt platforms and with different functionality

// LINUX is used to indicate that the code is compiled on a linux platform. Otherwise the platform is assumed to be windows 7
#define WINDOWS

// If DEBUG is defined then a number of files are created in the debugOut directory containing ASCII data sets
//#define DEBUG

// The data format for the hash_table - can be either int or double
#define HASH_TYPE double

// uncomment this line to increas the ammount of output to the screen
//#define OUTPUT

// Terminate on CUDA error?
#define CUDA_ERROR_TERMINATE

// Cuda device to use
#define CUDA_DEVICE 0