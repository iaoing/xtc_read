#ifndef PLFS_XTC_H
#define PLFS_XTC_H

#include <stdio.h>
#include <limits.h>

// Error codes for mdio_errno
#define MDIO_SUCCESS		0
#define MDIO_BADFORMAT		1
#define MDIO_EOF		2
#define MDIO_BADPARAMS		3
#define MDIO_IOERROR		4
#define MDIO_BADPRECISION	5
#define MDIO_BADMALLOC		6
#define MDIO_CANTOPEN		7
#define MDIO_BADEXTENSION	8
#define MDIO_UNKNOWNFMT		9
#define MDIO_CANTCLOSE		10
#define MDIO_WRONGFORMAT	11
#define MDIO_SIZEERROR		12
#define MDIO_MODEERROR		13
#define MDIO_UNKNOWNERROR	1000

#define XTC_READ	0
#define XTC_WRITE	1

#define MDIO_MAX_ERRVAL		11

static int mdio_errcode;	// Last error code

#define TRX_MAGIC	1993	// Magic number for .trX files
#define XTC_MAGIC	1995	// Magic number for .xtc files
#define MAX_GRO_LINE	500	// Maximum line length of .gro files
#define MAX_G96_LINE	500	// Maximum line length of .g96 files
#define MAX_TRX_TITLE	80	// Maximum length of a title in .trX
#define MAX_MDIO_TITLE	80	// Maximum supported title length
#define ANGS_PER_NM	10	// Unit conversion factor
#define ANGS2_PER_NM2	100	// Unit conversion factor

// All the supported file types and their respective extensions
#define MDFMT_GRO		1
#define MDFMT_TRR		2
#define MDFMT_G96		3
#define MDFMT_TRJ		4
#define MDFMT_XTC		5


#define BYTES_PER_XTC_UNIT  (4)
static char xtc_zero[BYTES_PER_XTC_UNIT] = {0, 0, 0, 0};

#define MAXABS INT_MAX-2

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// integer table used in decompression
static int xtc_magicints[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0,8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
	80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
	1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 16384,
	20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031, 131072,
	165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255,
	1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304,
	5284491, 6658042, 8388607, 10568983, 13316085, 16777216 };

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(xtc_magicints) / sizeof(*xtc_magicints))


void XTC_DEBUG(const char *file, const char *func, const int line) ;


// A generic i/o structure that contains information about the
// file itself and the input/output state
typedef struct {
	FILE *	f;	// Pointer to the file
	int	fmt;	// The file format
	int	prec;	// Real number precision
	int	rev;	// Reverse endiannism?
	int mode;	// open mode, 0 means XTC_READ, 1 means XTC_WRITE
} md_file;


// A format-independent structure to hold header data from files
typedef struct {
	char title[MAX_MDIO_TITLE + 1];
	int natoms;
	float timeval;
} md_header;

// A format-independent structure to hold unit cell data
typedef struct {
  float A, B, C, alpha, beta, gamma;
} md_box;

// Timestep information
typedef struct {
	float *pos;	// Position array (3 * natoms)
	//float *vel;	// Velocity array ** (VMD doesn't use this) **
	//float *f;	// Force array ** (VMD doesn't use this) **
	//float *box;	// Computational box ** (VMD doesn't use this) **
	int natoms;	// Number of atoms
	int step;	// Simulation step
	float time;	// Time of simulation
	float precision;
	md_box *box;
	float matri_box[3][3];
} md_ts;

// Error descriptions for mdio_errno
static const char *mdio_errdescs[] = {
	"no error",
	"file does not match format",
	"unexpected end-of-file reached",
	"function called with bad parameters",
	"file i/o error",
	"unsupported precision",
	"out of memory",
	"cannot open file",
	"bad file extension",
	"unknown file format",
	"cannot close file",
	"wrong file format for this function",
	"binary i/o error: sizeof(int) != 4",
	"unknown file open mode",
	NULL
};


// Error reporting functions
int mdio_errno(void);
const char *mdio_errmsg(int);
static int mdio_seterror(int);


// Open a molecular dynamics file. The second parameter specifies
// the format of the file. If it is zero, the format is determined
// from the file extension. the third argument (if given) decides
// whether to read (==0) or to write (!= 0).
// using a default argument set to read for backward compatibility.
md_file *xtc_open(const char *, const int fmt, const int rw);

// Closes a molecular dynamics file.
int xtc_close(md_file *);


// .xtc file functions
int xtc_read_header(md_file *, md_header *);
int xtc_read_frame(md_file *, md_ts *);
int xtc_write_frame(md_file *, md_ts *);

// .xtc file functions 
int xtc_int(md_file *, int *);
int xtc_float(md_file *, float *);
int xtc_data(md_file *, char *, int );

// .xtc file functions for read
int get_int(md_file *, int *);
int get_float(md_file *, float *);
int get_data(md_file *, char *, int );

// .xtc file functions for write
int put_int(md_file *, int *);
int put_float(md_file *, float *);
int put_data(md_file *, char *, int );

// some tools for put or get int/float/data
unsigned int xdr_ntohl(unsigned int x);
unsigned int xdr_htonl(unsigned int x);
unsigned int xdr_swapbytes(unsigned int x);

// Miscellaneous functions
int mdio_readbox(md_box *, float *, float *, float *);


// for 3dfcoord
// some tools
int xtc_sizeofint(int );
int xtc_sizeofints(int , unsigned int *);
// read tools
int xtc_receivebits(int *, int );
void xtc_receiveints(int *, const int , int , unsigned int *, int *);
// write tools
void xtc_sendbits(int *, int , int );
void xtc_sendints(int *, const int , const int , unsigned int *, unsigned int *);
// main function
int xtc_3dfcoord(md_file *, float *, int *, float *);
int xtc3dfcoord(md_file *, float *, int *, float *);

// other functions to get more info
int xtc_get_water_no(FILE *fp, int *water_index, int water_length, int *wn_index);

// write uncompressed frame info to a file;
int xtc_write_uc_frame(md_file *mf, md_ts *);
int xtc_uc_3dfcoord(md_file *mf, float *, int *, float *);
int xtc_uc_data(md_file *, char *, int );


#endif	

