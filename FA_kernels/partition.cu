#include <stdio.h>
#include <stdlib.h>
#include "partition.h"
//#include "conio.h"

/*  ----------------------  */
/*  ---- buildpackage ----  */
/*  ----------------------  */

int buildpackagei(int *grid, int GRDROWS, int GRDCOLS, int beginy, int beginx, int partrows, int partcols, int *package)
{
    int y, x;
    int endy, endx;
    int copytop, copybottom, copyleft, copyright;
    int begincopyat, endcopyat;
    int beginp, endp;
    int *next;
    
    endy = beginy + (partrows-1);
    endx = beginx + (partcols-1);
    

    next = package + 2*PKGCOLS + 2; /* jump to third element of third row, i.e. element at 2,2*/
    
    for (y = beginy; y <= endy; ++y)
    {
        for (x = beginx; x <= endx; ++x)
        {
            *next = grid[y*GRDCOLS + x];
            next++;
        }
        
        next += 4; /* to account for left and right edges */
    }
    
    copytop = copybottom = copyleft = copyright = 1;
    
    if (beginy == 0)         copytop    = 0;
    if (beginx == 0)         copyleft   = 0;
    if (endy   == GRDROWS-1) copybottom = 0;
    if (endx   == GRDCOLS-1) copyright  = 0;
    
    beginp = beginy*GRDCOLS + beginx;
    endp = endy*GRDCOLS + endx;
    
    if (copytop)
    {
        begincopyat = beginp - GRDCOLS                - copyleft*1;  /* move one row up |                      | may move one col to the left */
        endcopyat   = beginp - GRDCOLS + (partcols-1) + copyright*1; /* move one row up | jump to last element | may move one col to the right */
        next = package + PKGCOLS + 1 + (1-copyleft)*1; /* jump to first element of second row | may move two cols to the right */
        
        for (x = begincopyat; x <= endcopyat; ++x)
        {
            *next = grid[x];
            next++;
        }
        
        if (beginy >= 2)
        {
            begincopyat = beginp - 2*GRDCOLS                - copyleft*2;  /* move one row up |                      | may move one col to the left */
            endcopyat   = beginp - 2*GRDCOLS + (partcols-1) + copyright*2; /* move one row up | jump to last element | may move one col to the right */
            next = package + (1-copyleft)*2; /* jump to first element of first row | may move twp cols to the right */
            
            for (x = begincopyat; x <= endcopyat; ++x)
            {
                *next = grid[x];
                next++;
            }
        }
    }
    
    if (copybottom)
    {
        begincopyat = endp + GRDCOLS - (partcols-1) - copyleft*1;  /*  move one row down | jump to first element | may move one col to the left */
        endcopyat   = endp + GRDCOLS                + copyright*1; /*  move one row down |                       | may move one col to the right */
        next = package + (PKGROWS-2)*PKGCOLS + 1 + (1-copyleft)*1;    /* jump to first element of last row | may move one one to the right */
        
        for (x = begincopyat; x <= endcopyat; ++x)
        {
            *next = grid[x];
            next++;
        }
        
        if (endy <= GRDROWS-3)
        {
            begincopyat = endp + 2*GRDCOLS - (partcols-1) - copyleft*2;  /*  move one row down | jump to first element | may move one col to the left */
            endcopyat   = endp + 2*GRDCOLS                + copyright*2; /*  move one row down |                       | may move one col to the right */
            next = package + (PKGROWS-1)*PKGCOLS + (1-copyleft)*2;    /* jump to first element of last row | may move one one to the right */
            
            for (x = begincopyat; x <= endcopyat; ++x)
            {
                *next = grid[x];
                next++;
            }
        }
    }
    
    if (copyleft)
    {
        begincopyat = beginp - 1                - copytop*GRDCOLS;    /* move one col to the left |                       | may move one row up */
        endcopyat   =   endp - 1 - (partcols-1) + copybottom*GRDCOLS; /* move one col to the left | jump to first element | may move one row down */
        next = package + 1 + PKGCOLS + (1-copytop)*PKGCOLS;                        /* jump to first element of first row | may move one row down */
        
        for (x = begincopyat; x <= endcopyat; x += GRDCOLS)
        {
            *next = grid[x];
            next += PKGCOLS;
        }
        
        if (beginx >= 2)
        {
            begincopyat = beginp - 2                - copytop*2*GRDCOLS;    /* move one col to the left |                       | may move one row up */
            endcopyat   =   endp - 2 - (partcols-1) + copybottom*2*GRDCOLS; /* move one col to the left | jump to first element | may move one row down */
            next = package + (1-copytop)*2*PKGCOLS;                        /* jump to first element of first row | may move one row down */
            
            for (x = begincopyat; x <= endcopyat; x += GRDCOLS)
            {
                *next = grid[x];
                next += PKGCOLS;
            }
        }
    }
    
    if (copyright)
    {
        begincopyat = beginp + 1 + (partcols-1) - copytop*GRDCOLS;    /* move one col to the right | jump to last element | may move one row up */
        endcopyat   =   endp + 1                + copybottom*GRDCOLS; /* move one col to the right |                      | may move one row down */
        next = package - 1 + PKGCOLS + (PKGCOLS-1) + (1-copytop)*PKGCOLS;          /* jump to last element of first row | may move one row down */
        
        for (x = begincopyat; x <= endcopyat; x += GRDCOLS)
        {
            *next = grid[x];
            next += PKGCOLS;
        }
        
        if (endx <= GRDCOLS-3)
        {
            begincopyat = beginp + 2 + (partcols-1) - copytop*2*GRDCOLS;    /* move one col to the right | jump to last element | may move one row up */
            endcopyat   =   endp + 2                + copybottom*2*GRDCOLS; /* move one col to the right |                      | may move one row down */
            next = package + (PKGCOLS-1) + (1-copytop)*2*PKGCOLS;          /* jump to last element of first row | may move one row down */
            
            for (x = begincopyat; x <= endcopyat; x += GRDCOLS)
            {
                *next = grid[x];
                next += PKGCOLS;
            }
        }
    }
    
    return 0;
}

int buildpackagef(double *grid, int GRDROWS, int GRDCOLS, int beginy, int beginx, int partrows, int partcols, double *package)
{
    int y, x;
    int endy, endx;
    int copytop, copybottom, copyleft, copyright;
    int begincopyat, endcopyat;
    int beginp, endp;
    double *next;
    
    endy = beginy + (partrows-1);
    endx = beginx + (partcols-1);
    
    next = package + 2*PKGCOLS + 2; /* jump to third element of third row, i.e. element at 2,2*/
    

	//for (int s=0; s <10; s++ ){
	//	printf (" %lf \n", grid[s]);
	//}
	//_getch();




    for (y = beginy; y <= endy; ++y)
    {
        for (x = beginx; x <= endx; ++x)
        {
            *next = grid[y*GRDCOLS + x];
            next++;
        }
        
        next += 4; /* to account for left and right edges */
    }
    
    copytop = copybottom = copyleft = copyright = 1;
    
    if (beginy == 0)         copytop    = 0;
    if (beginx == 0)         copyleft   = 0;
    if (endy   == GRDROWS-1) copybottom = 0;
    if (endx   == GRDCOLS-1) copyright  = 0;
    
    beginp = beginy*GRDCOLS + beginx;
    endp = endy*GRDCOLS + endx;
    
    if (copytop)
    {
        begincopyat = beginp - GRDCOLS                - copyleft*1;  /* move one row up |                      | may move one col to the left */
        endcopyat   = beginp - GRDCOLS + (partcols-1) + copyright*1; /* move one row up | jump to last element | may move one col to the right */
        next = package + PKGCOLS + 1 + (1-copyleft)*1; /* jump to first element of second row | may move two cols to the right */
        
        for (x = begincopyat; x <= endcopyat; ++x)
        {
            *next = grid[x];
            next++;
        }
        
        if (beginy >= 2)
        {
            begincopyat = beginp - 2*GRDCOLS                - copyleft*2;  /* move one row up |                      | may move one col to the left */
            endcopyat   = beginp - 2*GRDCOLS + (partcols-1) + copyright*2; /* move one row up | jump to last element | may move one col to the right */
            next = package + (1-copyleft)*2; /* jump to first element of first row | may move twp cols to the right */
            
            for (x = begincopyat; x <= endcopyat; ++x)
            {
                *next = grid[x];
                next++;
            }
        }
    }
    
    if (copybottom)
    {
        begincopyat = endp + GRDCOLS - (partcols-1) - copyleft*1;  /*  move one row down | jump to first element | may move one col to the left */
        endcopyat   = endp + GRDCOLS                + copyright*1; /*  move one row down |                       | may move one col to the right */
        next = package + (PKGROWS-2)*PKGCOLS + 1 + (1-copyleft)*1;    /* jump to first element of last row | may move one one to the right */
        
        for (x = begincopyat; x <= endcopyat; ++x)
        {
            *next = grid[x];
            next++;
        }
        
        if (endy <= GRDROWS-3)
        {
            begincopyat = endp + 2*GRDCOLS - (partcols-1) - copyleft*2;  /*  move one row down | jump to first element | may move one col to the left */
            endcopyat   = endp + 2*GRDCOLS                + copyright*2; /*  move one row down |                       | may move one col to the right */
            next = package + (PKGROWS-1)*PKGCOLS + (1-copyleft)*2;    /* jump to first element of last row | may move one one to the right */
            
            for (x = begincopyat; x <= endcopyat; ++x)
            {
                *next = grid[x];
                next++;
            }
        }
    }
    
    if (copyleft)
    {
        begincopyat = beginp - 1                - copytop*GRDCOLS;    /* move one col to the left |                       | may move one row up */
        endcopyat   =   endp - 1 - (partcols-1) + copybottom*GRDCOLS; /* move one col to the left | jump to first element | may move one row down */
        next = package + 1 + PKGCOLS + (1-copytop)*PKGCOLS;                        /* jump to first element of first row | may move one row down */
        
        for (x = begincopyat; x <= endcopyat; x += GRDCOLS)
        {
            *next = grid[x];
            next += PKGCOLS;
        }
        
        if (beginx >= 2)
        {
            begincopyat = beginp - 2                - copytop*2*GRDCOLS;    /* move one col to the left |                       | may move one row up */
            endcopyat   =   endp - 2 - (partcols-1) + copybottom*2*GRDCOLS; /* move one col to the left | jump to first element | may move one row down */
            next = package + (1-copytop)*2*PKGCOLS;                        /* jump to first element of first row | may move one row down */
            
            for (x = begincopyat; x <= endcopyat; x += GRDCOLS)
            {
                *next = grid[x];
                next += PKGCOLS;
            }
        }
    }
    
    if (copyright)
    {
        begincopyat = beginp + 1 + (partcols-1) - copytop*GRDCOLS;    /* move one col to the right | jump to last element | may move one row up */
        endcopyat   =   endp + 1                + copybottom*GRDCOLS; /* move one col to the right |                      | may move one row down */
        next = package - 1 + PKGCOLS + (PKGCOLS-1) + (1-copytop)*PKGCOLS;          /* jump to last element of first row | may move one row down */
        
        for (x = begincopyat; x <= endcopyat; x += GRDCOLS)
        {
            *next = grid[x];
            next += PKGCOLS;
        }
        
        if (endx <= GRDCOLS-3)
        {
            begincopyat = beginp + 2 + (partcols-1) - copytop*2*GRDCOLS;    /* move one col to the right | jump to last element | may move one row up */
            endcopyat   =   endp + 2                + copybottom*2*GRDCOLS; /* move one col to the right |                      | may move one row down */
            next = package + (PKGCOLS-1) + (1-copytop)*2*PKGCOLS;          /* jump to last element of first row | may move one row down */
            
            for (x = begincopyat; x <= endcopyat; x += GRDCOLS)
            {
                *next = grid[x];
                next += PKGCOLS;
            }
        }
    }
    
    return 0;
}

/*  ---------------------------  */
/*  ---- unwrapandreinsert ----  */
/*  ---------------------------  */

int unwrapandreinserti(int *grid, int gridRows, int gridColumns, int beginy, int beginx, int partrows, int partcols, int *package)
{
    int y, x;
    int index;
    int gridy, gridx;
    
    gridy = beginy;
    gridx = beginx;

    for (y = 0; y < partrows; ++y)
    {
        for (x = 0; x < partcols; ++x)
        {
            index = y*(partcols+4) + x + 2*(partcols+4) + 2;
            /*if (index > (partcols +4) * (partrows + 4)) {
              printf("Invalid memory access!\n");
              exit(0);
            }*/
            /*if (gridy*gridColumns + gridx > gridRows * gridColumns) {
              printf("Invalid memory access2!\n");
              printf("%d > %d\n", gridy*gridColumns + gridx, gridRows * gridColumns);
              printf("gridy = %d, gridx = %d\n", gridy, gridx);
              printf("beginx = %d, partcols = %d\n", beginx, partcols);
              exit(0);
            }*/
            grid[gridy*gridColumns + gridx] = package[index];
            gridx++;
        }
        
        gridx = beginx;
        gridy++;
    }
    
    return 0;
}

int unwrapandreinsertf(double *grid, int gridRows, int gridColumns, int beginy, int beginx, int partrows, int partcols, double *package)
{
    int y, x;
    int index;
    int gridy, gridx;
	int check_count = 0;
    
    gridy = beginy;
    gridx = beginx;

    for (y = 0; y < partrows; ++y)
    {
        for (x = 0; x < partcols; ++x)
        {
            index = y*(partcols+4) + x + 2*(partcols+4) + 2;
            grid[gridy*gridColumns + gridx] = package[index];
			if (package[index] > 0) check_count ++;
            gridx++;
        }
        
        gridx = beginx;
        gridy++;
    }
    
	//printf("check count above zero for fa grid %d \n", check_count);

    return 0;
}

/*  ----------------------  */
/*  ---- cleanpackage ----  */
/*  ----------------------  */

int cleanpackagei(int *package, int partrows, int partcols)
{
    int y, x;
    
    for (x = 0; x < PKGCOLS; ++x)
    {
        package[x] = 0;
        package[(PKGROWS-1)*PKGCOLS + x] = 0;
    }
    
    for (x = 0; x < PKGCOLS; ++x)
    {
        package[PKGCOLS + x] = 0;
        package[(PKGROWS-2)*PKGCOLS + x] = 0;
    }
    
    for (y = 0; y < PKGROWS; ++y)
    {
        package[y*PKGCOLS] = 0;
        package[y*PKGCOLS + (PKGCOLS-1)] = 0;
    }
    
    for (y = 0; y < PKGROWS; ++y)
    {
        package[y*PKGCOLS + 1] = 0;
        package[y*PKGCOLS + (PKGCOLS-1) - 1] = 0;
    }
    
    return 0;
}

int cleanpackagef(double *package, int partrows, int partcols)
{
    int y, x;
    
    for (x = 0; x < PKGCOLS; ++x)
    {
        package[x] = 0;
        package[(PKGROWS-1)*PKGCOLS + x] = 0;
    }
    
    for (x = 0; x < PKGCOLS; ++x)
    {
        package[PKGCOLS + x] = 0;
        package[(PKGROWS-2)*PKGCOLS + x] = 0;
    }
    
    for (y = 0; y < PKGROWS; ++y)
    {
        package[y*PKGCOLS] = 0;
        package[y*PKGCOLS + (PKGCOLS-1)] = 0;
    }
    
    for (y = 0; y < PKGROWS; ++y)
    {
        package[y*PKGCOLS + 1] = 0;
        package[y*PKGCOLS + (PKGCOLS-1) - 1] = 0;
    }
    
    return 0;
}

/*  ---------------------------  */
/*  ---- testing functions ----  */
/*  ---------------------------  */

int alterpackage(int *package, int partrows, int partcols)
{
    int y, x;
    int index;
    
    for (y = 0; y < partrows; ++y)
    {
        for (x = 0; x < partcols; ++x)
        {
            index = y*PKGCOLS + x + 2*PKGCOLS + 2;
            package[index] = index; /* y*partcols + (x+1) */
        }
    }
    
    return 0;
}

int printpackage(int *package, int partrows, int partcols)
{
    int y, x;
    
    for (y = 0; y < PKGROWS; ++y)
    {
        for (x = 0; x < PKGCOLS; ++x)
        {
            printf("%4d", package[y*PKGCOLS + x]);
        }
        
        printf("\n");
    }
    
    return 0;
}
