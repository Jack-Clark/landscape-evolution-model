
#include "lem.h"
#include "Data.h"

#ifndef memory_devH
#define memory_devH


void setdevicespace_FD(Data* data, Data* device);
void setdevicespace_FA(Data* data, Data* device);
void setdevicespace_Process(Data* data, Data* device);

void cleardevicespace_FD(Data* data, Data* device);
void cleardevicespace_FA(Data* data, Data* device);
void cleardevicespace_Process(Data* data, Data* device);

int createDeviceSpace(Data* data, Data* device);
int clearDeviceSpace(Data* data, Data* device);

int copyMask(Data* data, Data* device);

int cleargrid(Data* data);

#endif
