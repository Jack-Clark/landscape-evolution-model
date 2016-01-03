#define PKGROWS (partrows+4)
#define PKGCOLS (partcols+4)
int buildpackagei(int *grid, int GRDROWS, int GRDCOLS, int beginy, int beginx, int partrows, int partcols, int *package);
int buildpackagef(double *grid, int GRDROWS, int GRDCOLS, int beginy, int beginx, int partrows, int partcols, double *package);
int unwrapandreinserti(int *grid, int gridRows, int gridColumns, int beginy, int beginx, int partrows, int partcols, int *package);
int unwrapandreinsertf(double *grid, int gridRows, int gridColumns, int beginy, int beginx, int partrows, int partcols, double *package);
int cleanpackagei(int *package, int partrows, int partcols);
int cleanpackagef(double *package, int partrows, int partcols);
int alterpackage(int *package, int partrows, int partcols);
int printpackage(int *package, int partrows, int partcols);