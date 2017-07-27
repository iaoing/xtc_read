#include "plfs_xtc.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

void XTC_DEBUG(const char *file, const char *func, const int line) 
{
	printf("File: %s, function: %s, line: %d\n", file, func, line);
}

// Returns the last error code reported by any of the mdio functions
int mdio_errno(void) {
	return mdio_errcode;
}


// Returns a textual message regarding an mdio error code
const char *mdio_errmsg(int n) {
	if (n < 0 || n > MDIO_MAX_ERRVAL) return (char *) "unknown error";
	else return mdio_errdescs[n];
}


// Sets the error code and returns an appropriate return value
// for the calling function to return to its parent
// return -1 means curred some wrong
int mdio_seterror(int code) {
	mdio_errcode = code;
	return code ? -1 : 0;
}


// Open a molecular dynamics file. The second parameter specifies
// the format of the file. If it is zero, the format is determined
// from the file extension.
md_file *xtc_open(const char *fn, const int fmt, const int rw) {
	md_file *mf;

	if (!fn) {
		mdio_seterror(MDIO_BADPARAMS);
		return NULL;
	}

	// Allocate memory
	mf = (md_file *) malloc(sizeof(md_file));
	if (!mf) {
		mdio_seterror(MDIO_BADMALLOC);
		return NULL;
	}

	// Finally, open the file
    if (rw)			// write mode
        mf->f = fopen(fn, "wb");
    else			// read mode
        mf->f = fopen(fn, "rb");


	// Check for opening error
	if (!mf->f) {
		free(mf);
		mdio_seterror(MDIO_CANTOPEN);
		return NULL;
	}

	// set file's type ---- xtc
	mf->fmt = MDFMT_XTC;
	// set file's open mode
	mf->mode = rw;

	// File is opened, we're all set!
	mdio_seterror(MDIO_SUCCESS);
	return mf;
}

// Closes a molecular dynamics file.
int xtc_close(md_file *mf) {
	if (!mf) return mdio_seterror(MDIO_BADPARAMS);
	if (fclose(mf->f) == EOF) return mdio_seterror(MDIO_CANTCLOSE);

	free(mf);

	return mdio_seterror(MDIO_SUCCESS);
}




// Reads the header of a file (format independent)
int xtc_read_header(md_file *mf, md_header *mdh) {
	int n;
	if (!mf || !mdh) return mdio_seterror(MDIO_BADPARAMS);
	if (!mf->f) return mdio_seterror(MDIO_BADPARAMS);

	switch (mf->fmt) {
		case MDFMT_XTC:
			// printf("before memset: %d\n", __LINE__);
			memset(mdh, 0, sizeof(md_header));
			// Check magic number
			if (xtc_int(mf, &n) < 0) return -1;
			if (n != XTC_MAGIC) return mdio_seterror(MDIO_BADFORMAT);

			// Get number of atoms
			if (xtc_int(mf, &n) < 0) return -1;
			mdh->natoms = n;
			rewind(mf->f);
			return 0;

		default:
			return mdio_seterror(MDIO_UNKNOWNFMT);
	}
}


// xtc_read_frame() - reads a timestep from an .xtc file.
int xtc_read_frame(md_file *mf, md_ts *ts) {
	mf->mode = XTC_READ;

	int n;
	float f, x[3], y[3], z[3];

	int size = 0; // explicitly initialized to zero.
	float precision;

	if (!mf || !ts) return mdio_seterror(MDIO_BADPARAMS);
	if (!mf->f) return mdio_seterror(MDIO_BADPARAMS);
	if (mf->fmt != MDFMT_XTC) return mdio_seterror(MDIO_WRONGFORMAT);

	// Check magic number
	if (xtc_int(mf, &n) < 0) return -1;
	if (n != XTC_MAGIC) return mdio_seterror(MDIO_BADFORMAT);

	// Get number of atoms
	if (xtc_int(mf, &n) < 0) return -1;
	ts->natoms = n;

	// Get the simulation step
	if (xtc_int(mf, &n) < 0) return -1;
	ts->step = n;

	// Get the time value
	if (xtc_float(mf, &f) < 0) return -1;
	ts->time = f;

	// Read the basis vectors of the box
  	if ( (xtc_float(mf, &x[0]) < 0) ||
         (xtc_float(mf, &x[1]) < 0) ||
         (xtc_float(mf, &x[2]) < 0) ||
         (xtc_float(mf, &y[0]) < 0) ||
         (xtc_float(mf, &y[1]) < 0) ||
         (xtc_float(mf, &y[2]) < 0) ||
         (xtc_float(mf, &z[0]) < 0) ||
         (xtc_float(mf, &z[1]) < 0) ||
         (xtc_float(mf, &z[2]) < 0) )
    	return -1;
    
    ts->matri_box[0][0] = x[0];
    ts->matri_box[0][1] = x[1];
    ts->matri_box[0][2] = x[2];
    ts->matri_box[1][0] = y[0];
    ts->matri_box[1][1] = y[1];
    ts->matri_box[1][2] = y[2];
    ts->matri_box[2][0] = z[0];
    ts->matri_box[2][1] = z[1];
    ts->matri_box[2][2] = z[2];
	// Allocate the box and convert the vectors.
	ts->box = (md_box *) malloc(sizeof(md_box));
	if (mdio_readbox(ts->box, x, y, z) < 0) {
		free(ts->box);
		ts->box = NULL;
		return -1;
	}

	ts->pos = (float *) malloc(sizeof(float) * 3 * ts->natoms);
	if (!ts->pos) return mdio_seterror(MDIO_BADMALLOC);
	n = xtc3dfcoord(mf, ts->pos, &size, &precision);
	ts->precision = precision;

	if (n < 0) return -1;

	/* Now we're left with the job of scaling... */
	for (n = 0; n < ts->natoms * 3; n++)
		ts->pos[n] *= ANGS_PER_NM;

	return mdio_seterror(MDIO_SUCCESS);
}

// xtc_write_frame() - writes a timestep to an .xtc file.
int xtc_write_frame(md_file *mf, md_ts *ts) {
	mf->mode = XTC_WRITE;

	int n;
	float f;

	int size = ts->natoms; // explicitly initialized to zero.
	float precision = ts->precision;

	if (!mf || !ts) return mdio_seterror(MDIO_BADPARAMS);
	if (!mf->f) return mdio_seterror(MDIO_BADPARAMS);
	if (mf->fmt != MDFMT_XTC) return mdio_seterror(MDIO_WRONGFORMAT);

	// write magic number
	n = XTC_MAGIC;
	if (xtc_int(mf, &n) < 0) return -1;

	// write number of atoms
	n = ts->natoms;
	if (xtc_int(mf, &n) < 0) return -1;

	// write the simulation step
	n = ts->step;
	if (xtc_int(mf, &n) < 0) return -1;

	// write the time value
	f = ts->time;
	if (xtc_float(mf, &f) < 0) return -1;

	// write the basis vectors of the box
  	if ( (xtc_float(mf, &ts->matri_box[0][0]) < 0) ||
         (xtc_float(mf, &ts->matri_box[0][1]) < 0) ||
         (xtc_float(mf, &ts->matri_box[0][2]) < 0) ||
         (xtc_float(mf, &ts->matri_box[1][0]) < 0) ||
         (xtc_float(mf, &ts->matri_box[1][1]) < 0) ||
         (xtc_float(mf, &ts->matri_box[1][2]) < 0) ||
         (xtc_float(mf, &ts->matri_box[2][0]) < 0) ||
         (xtc_float(mf, &ts->matri_box[2][1]) < 0) ||
         (xtc_float(mf, &ts->matri_box[2][2]) < 0) )
    	return -1;

    /* Now we're left with the job of scaling... */
	for (n = 0; n < ts->natoms * 3; n++)
		ts->pos[n] /= ANGS_PER_NM;

	if (!ts->pos) return mdio_seterror(MDIO_BADPARAMS);
	n = xtc3dfcoord(mf, ts->pos, &size, &precision);
	if (n < 0) return -1;

	return mdio_seterror(MDIO_SUCCESS);
}


// .xtc file functions
int xtc_int(md_file *mf, int *i){
	if(mf->mode == XTC_READ) 
		return get_int(mf, i);
	else if(mf->mode == XTC_WRITE)
		return put_int(mf, i);
	else
		return mdio_seterror(MDIO_MODEERROR);
}

int xtc_float(md_file *mf, float *f){
	if(mf->mode == XTC_READ) 
		return get_float(mf, f);
	else if(mf->mode == XTC_WRITE)
		return put_float(mf, f);
	else
		return mdio_seterror(MDIO_MODEERROR);
}

int xtc_data(md_file *mf, char *buf, int len){
	if(mf->mode == XTC_READ) 
		return get_data(mf, buf, len);
	else if(mf->mode == XTC_WRITE)
		return put_data(mf, buf, len);
	else
		return mdio_seterror(MDIO_MODEERROR);
}

// .xtc file functions for read
// get_int() - reads an integer from an xtc file
int get_int(md_file *mf, int *i) {
	int mycopy;
	unsigned char *c;
	c = (char*)&mycopy;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);
        // sanity check.
    if (sizeof(int) != 4) {
    	XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
    	return mdio_seterror(MDIO_SIZEERROR);
    }

	if (fread(c, 1, 4, mf->f) != 4) {
		if (feof(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_EOF);
		}
		else if (ferror(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_IOERROR);
		}
		else {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_UNKNOWNERROR);
		}
	}

	// if (i) *i = c[3] + (c[2] << 8) + (c[1] << 16) + (c[0] << 24);
	*i = (int)(xdr_ntohl (mycopy));
	return mdio_seterror(MDIO_SUCCESS);
}


// get_float() - reads a float from an xtc file
int get_float(md_file *mf, float *f) {
	unsigned char *c;
	int i, mycopy;
	c = (char*)&mycopy;

	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	if (fread(c, 1, 4, mf->f) != 4) {
		if (feof(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_EOF);
		}
		else if (ferror(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_IOERROR);
		}
		else {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_UNKNOWNERROR);
		}
	}

	if (f) {
		// By reading the number in as an integer and then
		// copying it to a floating point number we can
		// ensure proper endianness
		// i = c[3] + (c[2] << 8) + (c[1] << 16) + (c[0] << 24);
		i = xdr_ntohl (mycopy);
		memcpy(f, &i, 4);
	}
	return mdio_seterror(MDIO_SUCCESS);
}


// get_data() - reads a specific amount of data from an xtc
// file using the xdr format.
int get_data(md_file *mf, char *buf, int len) {
	if (!mf || len < 1) return mdio_seterror(MDIO_BADPARAMS);
	size_t slen = (size_t)len;
	if (buf) {
		if (fread(buf, 1, slen, mf->f) != slen) {
			if (feof(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_EOF);
			}
			else if (ferror(mf->f)) {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_IOERROR);
			}
			else {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_UNKNOWNERROR);
			}
		}
		if (len % 4) {
			if (fseek(mf->f, 4 - (len % 4), SEEK_CUR)) {
				if (feof(mf->f)) {
					XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
					return mdio_seterror(MDIO_EOF);
				}
				else if (ferror(mf->f)) {
					XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
					return mdio_seterror(MDIO_IOERROR);
				}
				else {
					XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
					return mdio_seterror(MDIO_UNKNOWNERROR);
				}
			}
		}
	}
	else {
		int newlen;
		newlen = len;
		if (len % 4) newlen += (4 - (len % 4));
		if (fseek(mf->f, newlen, SEEK_CUR)) {
			if (feof(mf->f)) {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_EOF);
			}
			else if (ferror(mf->f)) {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_IOERROR);
			}
			else {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_UNKNOWNERROR);
			}
		}
	}
	return len;
}

// .xtc file functions for write
// put_int() - writes an integer to an xtc file
int put_int(md_file *mf, int *i) {
	if (!mf) return mdio_seterror(MDIO_BADPARAMS);
        // sanity check.
    if (sizeof(int) != 4) return mdio_seterror(MDIO_SIZEERROR);

    int mycopy = xdr_htonl (*i);
 
	if (fwrite((char *)&mycopy, 1, 4, mf->f) != 4) {
		if (feof(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_EOF);
		}
		else if (ferror(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_IOERROR);
		}
		else {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_UNKNOWNERROR);
		}
	}
	return mdio_seterror(MDIO_SUCCESS);
}


// put_float() - writes a float to an xtc file
int put_float(md_file *mf, float *f) {
	if (!mf) return mdio_seterror(MDIO_BADPARAMS);

	int mycopy = xdr_htonl (*(int*)f);
	
	if (fwrite((char *)&mycopy, 4, 1, mf->f) != 1) {
		if (feof(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_EOF);
		}
		else if (ferror(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_IOERROR);
		}
		else {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_UNKNOWNERROR);
		}
	}
	return mdio_seterror(MDIO_SUCCESS);
}


// put_data() - writes a specific amount of data to an xtc
// file using the xdr format.
int put_data(md_file *mf, char *buf, int len) {
	if (!mf || !buf || len < 1) return mdio_seterror(MDIO_BADPARAMS);
	size_t slen = (size_t)len;

	if (fwrite(buf, 1, slen, mf->f) != slen) {
		if (feof(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_EOF);
		}
		else if (ferror(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_IOERROR);
		}
		else {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_UNKNOWNERROR);
		}
	}

	if (len % 4) {
		// round bytes count to full xtc/xdr units
		int rndup = BYTES_PER_XTC_UNIT - (len % 4);
		if (fwrite(xtc_zero, rndup, 1, mf->f) != 1) {
			if (feof(mf->f)) {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_EOF);
			}
			else if (ferror(mf->f)) {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_IOERROR);
			}
			else {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_UNKNOWNERROR);
			}
		}
	}
	return len;
}

// some tools for get for put data
// xdr_ntohl -> (*char)int to real int for read
unsigned int xdr_ntohl(unsigned int x)
{
    short s = 0x0F00;
    if (*((char *)&s) == (char)0x0F)
    {
        /* bigendian, do nothing */
        return x;
    }
    else
    {
        /* smallendian, swap bytes */
        return xdr_swapbytes(x);
    }
}

// xdr_htonl -> real int to (char*)int for write
unsigned int xdr_htonl(unsigned int x)
{
    short s = 0x0F00;
    if (*((char *)&s) == (char)0x0F)
    {
        /* bigendian, do nothing */
        return x;
    }
    else
    {
        /* smallendian,swap bytes */
        return xdr_swapbytes(x);
    }
}

unsigned int xdr_swapbytes(unsigned int x)
{
    unsigned int y;
    int i;
    char        *px = (char *)&x;
    char        *py = (char *)&y;

    for (i = 0; i < 4; i++)
    {
        py[i] = px[3-i];
    }

    return y;
}




// Converts box basis vectors to A, B, C, alpha, beta, and gamma.  
// Stores values in md_box struct, which should be allocated before calling
// this function.
int mdio_readbox(md_box *box, float *x, float *y, float *z) {
  float A, B, C;

  if (!box) {
    return mdio_seterror(MDIO_BADPARAMS);
  }

  // A, B, C are the lengths of the x, y, z vectors, respectively
  A = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ) * ANGS_PER_NM;
  B = sqrt( y[0]*y[0] + y[1]*y[1] + y[2]*y[2] ) * ANGS_PER_NM;
  C = sqrt( z[0]*z[0] + z[1]*z[1] + z[2]*z[2] ) * ANGS_PER_NM;
  if ((A<=0) || (B<=0) || (C<=0)) {
    /* Use zero-length box size and set angles to 90. */
    box->A = box->B = box->C = 0;
    box->alpha = box->beta = box->gamma = 90;
  } else {
    box->A = A;
    box->B = B;
    box->C = C;
  
    // gamma, beta, alpha are the angles between the x & y, x & z, y & z
    // vectors, respectively
    box->gamma = acos( (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])*ANGS2_PER_NM2/(A*B) ) * 90.0/M_PI_2;
    box->beta = acos( (x[0]*z[0]+x[1]*z[1]+x[2]*z[2])*ANGS2_PER_NM2/(A*C) ) * 90.0/M_PI_2;
    box->alpha = acos( (y[0]*z[0]+y[1]*z[1]+y[2]*z[2])*ANGS2_PER_NM2/(B*C) ) * 90.0/M_PI_2; 
  }
  return mdio_seterror(MDIO_SUCCESS);
}



///////////////////////////////////////////////////////////////////////
// This algorithm is an implementation of the 3dfcoord algorithm
// written by Frans van Hoesel (hoesel@chem.rug.nl) as part of the
// Europort project in 1995.
///////////////////////////////////////////////////////////////////////

// returns the number of bits in the binary expansion of
// the given integer.
int xtc_sizeofint(int size) {
	unsigned int num = 1;
	unsigned int ssize = (unsigned int)size;
	int nbits = 0;

	while (ssize >= num && nbits < 32) {
		nbits++;
		num <<= 1;
	}
	return nbits;
}

// calculates the number of bits a set of integers, when compressed,
// will take up.
int xtc_sizeofints(int nints, unsigned int *sizes) {
	int i;
	unsigned int num;
	unsigned int nbytes, nbits, bytes[32], bytecnt, tmp;
	nbytes = 1;
	bytes[0] = 1;
	nbits = 0;
	for (i=0; i < nints; i++) {	
		tmp = 0;
		for (bytecnt = 0; bytecnt < nbytes; bytecnt++) {
			tmp = bytes[bytecnt] * sizes[i] + tmp;
			bytes[bytecnt] = tmp & 0xff;
			tmp >>= 8;
		}
		while (tmp != 0) {
			bytes[bytecnt++] = tmp & 0xff;
			tmp >>= 8;
		}
		nbytes = bytecnt;
	}
	num = 1;
	nbytes--;
	while (bytes[nbytes] >= num) {
		nbits++;
		num *= 2;
	}
	return nbits + nbytes * 8;
}

// reads bits from a buffer.    
int xtc_receivebits(int *buf, int nbits) {
	int cnt, num; 
	unsigned int lastbits, lastbyte;
	unsigned char * cbuf;
	int mask = (1 << nbits) -1;

	cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
	cnt = buf[0];
	lastbits = (unsigned int) buf[1];
	lastbyte = (unsigned int) buf[2];

	num = 0;
	while (nbits >= 8) {
		lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
		num |=  (lastbyte >> lastbits) << (nbits - 8);
		nbits -=8;
	}
	if (nbits > 0) {
		if (lastbits < (unsigned int)nbits) {
			lastbits += 8;
			lastbyte = (lastbyte << 8) | cbuf[cnt++];
		}
		lastbits -= nbits;
		num |= (lastbyte >> lastbits) & ((1 << nbits) -1);
	}
	num &= mask;
	buf[0] = cnt;
	buf[1] = lastbits;
	buf[2] = lastbyte;
	return num; 
}

// decompresses small integers from the buffer
// sizes parameter has to be non-zero to prevent divide-by-zero
void xtc_receiveints(int *buf, const int nints, int nbits,
			unsigned int *sizes, int *nums) {
	int bytes[32];
	int i, j, nbytes, p, num;

	bytes[1] = bytes[2] = bytes[3] = 0;
	nbytes = 0;
	while (nbits > 8) {
		bytes[nbytes++] = xtc_receivebits(buf, 8);
		nbits -= 8;
	}
	if (nbits > 0) {
		bytes[nbytes++] = xtc_receivebits(buf, nbits);
	}
	for (i = nints-1; i > 0; i--) {
		num = 0;
		for (j = nbytes-1; j >=0; j--) {
			num = (num << 8) | bytes[j];
			p = num / sizes[i];
			bytes[j] = p;
			num = num - p * sizes[i];
		}
		nums[i] = num;
	}
	nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}


/*____________________________________________________________________________
 |
 | sendbits - encode num into buf using the specified number of bits
 |
 | This routines appends the value of num to the bits already present in
 | the array buf. You need to give it the number of bits to use and you
 | better make sure that this number of bits is enough to hold the value
 | Also num must be positive.
 |
 */

void xtc_sendbits(int buf[], int num_of_bits, int num)
{

    unsigned int    cnt, lastbyte;
    int             lastbits;
    unsigned char * cbuf;

    cbuf     = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt      = (unsigned int) buf[0];
    lastbits = buf[1];
    lastbyte = (unsigned int) buf[2];
    while (num_of_bits >= 8)
    {
        lastbyte     = (lastbyte << 8) | ((num >> (num_of_bits -8)) /* & 0xff*/);
        cbuf[cnt++]  = lastbyte >> lastbits;
        num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
        lastbyte  = (lastbyte << num_of_bits) | num;
        lastbits += num_of_bits;
        if (lastbits >= 8)
        {
            lastbits   -= 8;
            cbuf[cnt++] = lastbyte >> lastbits;
        }
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits > 0)
    {
        cbuf[cnt] = lastbyte << (8 - lastbits);
    }
}


/*____________________________________________________________________________
 |
 | sendints - send a small set of small integers in compressed format
 |
 | this routine is used internally by xdr3dfcoord, to send a set of
 | small integers to the buffer.
 | Multiplication with fixed (specified maximum ) sizes is used to get
 | to one big, multibyte integer. Allthough the routine could be
 | modified to handle sizes bigger than 16777216, or more than just
 | a few integers, this is not done, because the gain in compression
 | isn't worth the effort. Note that overflowing the multiplication
 | or the byte buffer (32 bytes) is unchecked and causes bad results.
 |
 */

void xtc_sendints(int buf[], const int num_of_ints, const int num_of_bits,
                     unsigned int sizes[], unsigned int nums[])
{

    int          i, num_of_bytes, bytecnt;
    unsigned int bytes[32], tmp;

    tmp          = nums[0];
    num_of_bytes = 0;
    do
    {
        bytes[num_of_bytes++] = tmp & 0xff;
        tmp                 >>= 8;
    }
    while (tmp != 0);

    for (i = 1; i < num_of_ints; i++)
    {
        if (nums[i] >= sizes[i])
        {
            fprintf(stderr, "major breakdown in sendints num %u doesn't "
                    "match size %u\n", nums[i], sizes[i]);
            exit(1);
        }
        /* use one step multiply */
        tmp = nums[i];
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
            tmp            = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp          >>= 8;
        }
        while (tmp != 0)
        {
            bytes[bytecnt++] = tmp & 0xff;
            tmp            >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8)
    {
        for (i = 0; i < num_of_bytes; i++)
        {
            xtc_sendbits(buf, 8, bytes[i]);
        }
        xtc_sendbits(buf, num_of_bits - num_of_bytes * 8, 0);
    }
    else
    {
        for (i = 0; i < num_of_bytes-1; i++)
        {
            xtc_sendbits(buf, 8, bytes[i]);
        }
        xtc_sendbits(buf, num_of_bits- (num_of_bytes -1) * 8, bytes[i]);
    }
}



// function that actually reads and writes compressed coordinates  
// only sport read feature
int xtc_3dfcoord(md_file *mf, float *fp, int *size, float *precision) {
	static int *ip = NULL;
	static int oldsize;
	static int *buf;

	int minint[3], maxint[3], *lip;
	int smallidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
	int flag, k;
	int small, smaller, i, is_smaller, run;
	float *lfp;
	int tmp, *thiscoord,  prevcoord[3];

	int bufsize, lsize;
	unsigned int bitsize;
	float inv_precision;

        /* avoid uninitialized data compiler warnings */
        bitsizeint[0] = 0;
        bitsizeint[1] = 0;
        bitsizeint[2] = 0;

	if (xtc_int(mf, &lsize) < 0) return -1;

	if (*size != 0 && lsize != *size) return mdio_seterror(MDIO_BADFORMAT);

	*size = lsize;
	size3 = *size * 3;
	if (*size <= 9) {
		for (i = 0; i < *size; i++) {
			if (xtc_float(mf, fp + (3 * i)) < 0) return -1;
			if (xtc_float(mf, fp + (3 * i) + 1) < 0) return -1;
			if (xtc_float(mf, fp + (3 * i) + 2) < 0) return -1;
		}
		return *size;
	}
	xtc_float(mf, precision);
	if (ip == NULL) {
		ip = (int *)malloc(size3 * sizeof(*ip));
		if (ip == NULL) return mdio_seterror(MDIO_BADMALLOC);
		bufsize = (int) (size3 * 1.2);
		buf = (int *)malloc(bufsize * sizeof(*buf));
		if (buf == NULL) return mdio_seterror(MDIO_BADMALLOC);
		oldsize = *size;
	} else if (*size > oldsize) {
		ip = (int *)realloc(ip, size3 * sizeof(*ip));
		if (ip == NULL) return mdio_seterror(MDIO_BADMALLOC);
		bufsize = (int) (size3 * 1.2);
		buf = (int *)realloc(buf, bufsize * sizeof(*buf));
		if (buf == NULL) return mdio_seterror(MDIO_BADMALLOC);
		oldsize = *size;
	}
	buf[0] = buf[1] = buf[2] = 0;

	xtc_int(mf, &(minint[0]));
	xtc_int(mf, &(minint[1]));
	xtc_int(mf, &(minint[2]));

	xtc_int(mf, &(maxint[0]));
	xtc_int(mf, &(maxint[1]));
	xtc_int(mf, &(maxint[2]));
		
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
		bitsizeint[0] = xtc_sizeofint(sizeint[0]);
		bitsizeint[1] = xtc_sizeofint(sizeint[1]);
		bitsizeint[2] = xtc_sizeofint(sizeint[2]);
		bitsize = 0; /* flag the use of large sizes */
	} else {
		bitsize = xtc_sizeofints(3, sizeint);
	}

	xtc_int(mf, &smallidx);
	smaller = xtc_magicints[FIRSTIDX > smallidx - 1 ? FIRSTIDX : smallidx - 1] / 2;
	small = xtc_magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx];

	/* check for zero values that would yield corrupted data */
	if ( !sizesmall[0] || !sizesmall[1] || !sizesmall[2] ) {
		printf("XTC corrupted, sizesmall==0 (case 1)\n");
		return -1;
	}


	/* buf[0] holds the length in bytes */
	if (xtc_int(mf, &(buf[0])) < 0) return -1;

	if (xtc_data(mf, (char *) &buf[3], (int) buf[0]) < 0) return -1;

	buf[0] = buf[1] = buf[2] = 0;

	lfp = fp;
	inv_precision = 1.0f / (*precision);
	run = 0;
	i = 0;
	lip = ip;
	while (i < lsize) {
		thiscoord = (int *)(lip) + i * 3;

		if (bitsize == 0) {
			thiscoord[0] = xtc_receivebits(buf, bitsizeint[0]);
			thiscoord[1] = xtc_receivebits(buf, bitsizeint[1]);
			thiscoord[2] = xtc_receivebits(buf, bitsizeint[2]);
		} else {
			xtc_receiveints(buf, 3, bitsize, sizeint, thiscoord);
		}

		i++;
		thiscoord[0] += minint[0];
		thiscoord[1] += minint[1];
		thiscoord[2] += minint[2];

		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];
 

		flag = xtc_receivebits(buf, 1);
		is_smaller = 0;
		if (flag == 1) {
			run = xtc_receivebits(buf, 5);
			is_smaller = run % 3;
			run -= is_smaller;
			is_smaller--;
		}
		if (run > 0) {
			thiscoord += 3;
			for (k = 0; k < run; k+=3) {
				xtc_receiveints(buf, 3, smallidx, sizesmall, thiscoord);
				i++;
				thiscoord[0] += prevcoord[0] - small;
				thiscoord[1] += prevcoord[1] - small;
				thiscoord[2] += prevcoord[2] - small;
				if (k == 0) {
					/* interchange first with second atom for better
					 * compression of water molecules
					 */
					tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
					prevcoord[0] = tmp;
					tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
					prevcoord[1] = tmp;
					tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
					prevcoord[2] = tmp;
					*lfp++ = prevcoord[0] * inv_precision;
					*lfp++ = prevcoord[1] * inv_precision;
					*lfp++ = prevcoord[2] * inv_precision;

					if ( !sizesmall[0] || !sizesmall[1] || !sizesmall[2] ) {
						printf("XTC corrupted, sizesmall==0 (case 2)\n");
						return -1;
					}

				} else {
					prevcoord[0] = thiscoord[0];
					prevcoord[1] = thiscoord[1];
					prevcoord[2] = thiscoord[2];
				}
				*lfp++ = thiscoord[0] * inv_precision;
				*lfp++ = thiscoord[1] * inv_precision;
				*lfp++ = thiscoord[2] * inv_precision;
			}
		} else {
			*lfp++ = thiscoord[0] * inv_precision;
			*lfp++ = thiscoord[1] * inv_precision;
			*lfp++ = thiscoord[2] * inv_precision;		
		}
		smallidx += is_smaller;
		if (is_smaller < 0) {
			small = smaller;
			if (smallidx > FIRSTIDX) {
				smaller = xtc_magicints[smallidx - 1] /2;
			} else {
				smaller = 0;
			}
		} else if (is_smaller > 0) {
			smaller = small;
			small = xtc_magicints[smallidx] / 2;
		}
		sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx] ;
	}
	return 1;
}

/*____________________________________________________________________________
 |
 | xdr3dfcoord - read or write compressed 3d coordinates to xdr file.
 |
 | this routine reads or writes (depending on how you opened the file with
 | xdropen() ) a large number of 3d coordinates (stored in *fp).
 | The number of coordinates triplets to write is given by *size. On
 | read this number may be zero, in which case it reads as many as were written
 | or it may specify the number if triplets to read (which should match the
 | number written).
 | Compression is achieved by first converting all floating numbers to integer
 | using multiplication by *precision and rounding to the nearest integer.
 | Then the minimum and maximum value are calculated to determine the range.
 | The limited range of integers so found, is used to compress the coordinates.
 | In addition the differences between succesive coordinates is calculated.
 | If the difference happens to be 'small' then only the difference is saved,
 | compressing the data even more. The notion of 'small' is changed dynamically
 | and is enlarged or reduced whenever needed or possible.
 | Extra compression is achieved in the case of GROMOS and coordinates of
 | water molecules. GROMOS first writes out the Oxygen position, followed by
 | the two hydrogens. In order to make the differences smaller (and thereby
 | compression the data better) the order is changed into first one hydrogen
 | then the oxygen, followed by the other hydrogen. This is rather special, but
 | it shouldn't harm in the general case.
 |
 */

int xtc3dfcoord(md_file *mf, float *fp, int *size, float *precision)
{
    int     *ip  = NULL;
    int     *buf = NULL;
    int 	 bRead;

    /* preallocate a small buffer and ip on the stack - if we need more
       we can always malloc(). This is faster for small values of size: */
    unsigned     prealloc_size = 3*16;
    int          prealloc_ip[3*16], prealloc_buf[3*20];
    int          we_should_free = 0;

    int          minint[3], maxint[3], mindiff, *lip, diff;
    int          lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int          minidx, maxidx;
    unsigned     sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
    int          flag, k;
    int          smallnum, smaller, larger, i, is_small, is_smaller, run, prevrun;
    float       *lfp, lf;
    int          tmp, *thiscoord,  prevcoord[3];
    unsigned int tmpcoord[30];

    int          bufsize, xdrid, lsize;
    unsigned int bitsize;
    float        inv_precision;
    int          errval = 1;
    int          rc;

    bRead         = (mf->mode == XTC_READ);
    bitsizeint[0] = bitsizeint[1] = bitsizeint[2] = 0;
    prevcoord[0]  = prevcoord[1]  = prevcoord[2]  = 0;

    if (!bRead)
    {
        /* mf is open for writing */
    	// write atoms
        if(xtc_int(mf, size) < 0)
        {
        	return -1;
        }
        size3 = *size * 3;
        /* when the number of coordinates is small, don't try to compress; just
         * write them as floats using xdr_vector
         */
        if (*size <= 9)
        {
            for (i = 0; i < *size; i++) {
				if (xtc_float(mf, fp + (3 * i)) < 0) return -1;
				if (xtc_float(mf, fp + (3 * i) + 1) < 0) return -1;
				if (xtc_float(mf, fp + (3 * i) + 2) < 0) return -1;
			}
			return mdio_seterror(MDIO_SUCCESS);
        }

        if (xtc_float(mf, precision) < 0)
        {
            return -1;
        }

        if (size3 <= prealloc_size)
        {
            ip  = prealloc_ip;
            buf = prealloc_buf;
        }
        else
        {
            we_should_free = 1;
            bufsize        = size3 * 1.2;
            ip             = (int *)malloc((size_t)(size3 * sizeof(*ip)));
            buf            = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
            if (ip == NULL || buf == NULL)
            {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }
        /* buf[0-2] are special and do not contain actual data */
        buf[0]    = buf[1] = buf[2] = 0;
        minint[0] = minint[1] = minint[2] = INT_MAX;
        maxint[0] = maxint[1] = maxint[2] = INT_MIN;
        prevrun   = -1;
        lfp       = fp;
        lip       = ip;
        mindiff   = INT_MAX;
        oldlint1  = oldlint2 = oldlint3 = 0;
        while (lfp < fp + size3)
        {
            /* find nearest integer */
            if (*lfp >= 0.0)
            {
                lf = *lfp * *precision + 0.5;
            }
            else
            {
                lf = *lfp * *precision - 0.5;
            }
            if (fabs(lf) > MAXABS)
            {
                /* scaling would cause overflow */
                errval = 0;
            }
            lint1 = lf;
            if (lint1 < minint[0])
            {
                minint[0] = lint1;
            }
            if (lint1 > maxint[0])
            {
                maxint[0] = lint1;
            }
            *lip++ = lint1;
            lfp++;
            if (*lfp >= 0.0)
            {
                lf = *lfp * *precision + 0.5;
            }
            else
            {
                lf = *lfp * *precision - 0.5;
            }
            if (fabs(lf) > MAXABS)
            {
                /* scaling would cause overflow */
                errval = 0;
            }
            lint2 = lf;
            if (lint2 < minint[1])
            {
                minint[1] = lint2;
            }
            if (lint2 > maxint[1])
            {
                maxint[1] = lint2;
            }
            *lip++ = lint2;
            lfp++;
            if (*lfp >= 0.0)
            {
                lf = *lfp * *precision + 0.5;
            }
            else
            {
                lf = *lfp * *precision - 0.5;
            }
            if (fabs(lf) > MAXABS)
            {
                /* scaling would cause overflow */
                errval = 0;
            }
            lint3 = lf;
            if (lint3 < minint[2])
            {
                minint[2] = lint3;
            }
            if (lint3 > maxint[2])
            {
                maxint[2] = lint3;
            }
            *lip++ = lint3;
            lfp++;
            diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
            if (diff < mindiff && lfp > fp + 3)
            {
                mindiff = diff;
            }
            oldlint1 = lint1;
            oldlint2 = lint2;
            oldlint3 = lint3;
        }
        if ( (xtc_int(mf, &(minint[0])) < 0) ||
             (xtc_int(mf, &(minint[1])) < 0) ||
             (xtc_int(mf, &(minint[2])) < 0) ||
             (xtc_int(mf, &(maxint[0])) < 0) ||
             (xtc_int(mf, &(maxint[1])) < 0) ||
             (xtc_int(mf, &(maxint[2])) < 0))
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return -1;
        }

        if ((float)maxint[0] - (float)minint[0] >= MAXABS ||
            (float)maxint[1] - (float)minint[1] >= MAXABS ||
            (float)maxint[2] - (float)minint[2] >= MAXABS)
        {
            /* turning value in unsigned by subtracting minint
             * would cause overflow
             */
            errval = 0;
        }
        sizeint[0] = maxint[0] - minint[0]+1;
        sizeint[1] = maxint[1] - minint[1]+1;
        sizeint[2] = maxint[2] - minint[2]+1;

        /* check if one of the sizes is to big to be multiplied */
        if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
        {
            bitsizeint[0] = xtc_sizeofint(sizeint[0]);
            bitsizeint[1] = xtc_sizeofint(sizeint[1]);
            bitsizeint[2] = xtc_sizeofint(sizeint[2]);
            bitsize       = 0; /* flag the use of large sizes */
        }
        else
        {
            bitsize = xtc_sizeofints(3, sizeint);
        }
        lip      = ip;
        luip     = (unsigned int *) ip;
        smallidx = FIRSTIDX;
        while (smallidx < LASTIDX && xtc_magicints[smallidx] < mindiff)
        {
            smallidx++;
        }
        if (xtc_int(mf, &smallidx) < 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return -1;
        }

        maxidx       = MIN(LASTIDX, smallidx + 8);
        minidx       = maxidx - 8; /* often this equal smallidx */
        smaller      = xtc_magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
        smallnum     = xtc_magicints[smallidx] / 2;
        sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx];
        larger       = xtc_magicints[maxidx] / 2;
        i            = 0;
        while (i < *size)
        {
            is_small  = 0;
            thiscoord = (int *)(luip) + i * 3;
            if (smallidx < maxidx && i >= 1 &&
                abs(thiscoord[0] - prevcoord[0]) < larger &&
                abs(thiscoord[1] - prevcoord[1]) < larger &&
                abs(thiscoord[2] - prevcoord[2]) < larger)
            {
                is_smaller = 1;
            }
            else if (smallidx > minidx)
            {
                is_smaller = -1;
            }
            else
            {
                is_smaller = 0;
            }
            if (i + 1 < *size)
            {
                if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
                    abs(thiscoord[1] - thiscoord[4]) < smallnum &&
                    abs(thiscoord[2] - thiscoord[5]) < smallnum)
                {
                    /* interchange first with second atom for better
                     * compression of water molecules
                     */
                    tmp          = thiscoord[0]; thiscoord[0] = thiscoord[3];
                    thiscoord[3] = tmp;
                    tmp          = thiscoord[1]; thiscoord[1] = thiscoord[4];
                    thiscoord[4] = tmp;
                    tmp          = thiscoord[2]; thiscoord[2] = thiscoord[5];
                    thiscoord[5] = tmp;
                    is_small     = 1;
                }

            }
            tmpcoord[0] = thiscoord[0] - minint[0];
            tmpcoord[1] = thiscoord[1] - minint[1];
            tmpcoord[2] = thiscoord[2] - minint[2];
            if (bitsize == 0)
            {
                xtc_sendbits(buf, bitsizeint[0], tmpcoord[0]);
                xtc_sendbits(buf, bitsizeint[1], tmpcoord[1]);
                xtc_sendbits(buf, bitsizeint[2], tmpcoord[2]);
            }
            else
            {
                xtc_sendints(buf, 3, bitsize, sizeint, tmpcoord);
            }
            prevcoord[0] = thiscoord[0];
            prevcoord[1] = thiscoord[1];
            prevcoord[2] = thiscoord[2];
            thiscoord    = thiscoord + 3;
            i++;

            run = 0;
            if (is_small == 0 && is_smaller == -1)
            {
                is_smaller = 0;
            }
            while (is_small && run < 8*3)
            {
                if (is_smaller == -1 && (
                        SQR(thiscoord[0] - prevcoord[0]) +
                        SQR(thiscoord[1] - prevcoord[1]) +
                        SQR(thiscoord[2] - prevcoord[2]) >= smaller * smaller))
                {
                    is_smaller = 0;
                }

                tmpcoord[run++] = thiscoord[0] - prevcoord[0] + smallnum;
                tmpcoord[run++] = thiscoord[1] - prevcoord[1] + smallnum;
                tmpcoord[run++] = thiscoord[2] - prevcoord[2] + smallnum;

                prevcoord[0] = thiscoord[0];
                prevcoord[1] = thiscoord[1];
                prevcoord[2] = thiscoord[2];

                i++;
                thiscoord = thiscoord + 3;
                is_small  = 0;
                if (i < *size &&
                    abs(thiscoord[0] - prevcoord[0]) < smallnum &&
                    abs(thiscoord[1] - prevcoord[1]) < smallnum &&
                    abs(thiscoord[2] - prevcoord[2]) < smallnum)
                {
                    is_small = 1;
                }
            }
            if (run != prevrun || is_smaller != 0)
            {
                prevrun = run;
                xtc_sendbits(buf, 1, 1); /* flag the change in run-length */
                xtc_sendbits(buf, 5, run+is_smaller+1);
            }
            else
            {
                xtc_sendbits(buf, 1, 0); /* flag the fact that runlength did not change */
            }
            for (k = 0; k < run; k += 3)
            {
                xtc_sendints(buf, 3, smallidx, sizesmall, &tmpcoord[k]);
            }
            if (is_smaller != 0)
            {
                smallidx += is_smaller;
                if (is_smaller < 0)
                {
                    smallnum = smaller;
                    smaller  = xtc_magicints[smallidx-1] / 2;
                }
                else
                {
                    smaller  = smallnum;
                    smallnum = xtc_magicints[smallidx] / 2;
                }
                sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx];
            }
        }
        if (buf[1] != 0)
        {
            buf[0]++;
        }
        /* buf[0] holds the length in bytes */
        if (xtc_int(mf, &(buf[0])) < 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return -1;
        }


        rc = errval * (xtc_data(mf, (char *)&(buf[3]), (int)buf[0]));
        if (we_should_free)
        {
            free(ip);
            free(buf);
        }
        return mdio_seterror(MDIO_SUCCESS);

    }
    else
    {

        /* mf is open for reading */
        // return xtc_3dfcoord(mf, fp, size, precision);

    	// read natoms
        if (xtc_int(mf, &lsize) < 0)
        {
            return -1;
        }
        if (*size != 0 && lsize != *size)
        {
            fprintf(stderr, "wrong number of coordinates in xdr3dfcoord; "
                    "%d arg vs %d in file", *size, lsize);
            return mdio_seterror(MDIO_BADFORMAT);
        }
        *size = lsize;
        size3 = *size * 3;
        if (*size <= 9)
        {
            for (i = 0; i < *size; i++) {
				if (xtc_float(mf, fp + (3 * i)) < 0) return -1;
				if (xtc_float(mf, fp + (3 * i) + 1) < 0) return -1;
				if (xtc_float(mf, fp + (3 * i) + 2) < 0) return -1;
			}
			return *size;
        }
        if(xtc_float(mf, precision) < 0)
        {
        	return -1;
        }

        if (size3 <= prealloc_size)
        {
            ip  = prealloc_ip;
            buf = prealloc_buf;
        }
        else
        {
            we_should_free = 1;
            bufsize        = size3 * 1.2;
            ip             = (int *)malloc((size_t)(size3 * sizeof(*ip)));
            buf            = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
            if (ip == NULL || buf == NULL)
            {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }

        buf[0] = buf[1] = buf[2] = 0;

        if ( (xtc_int(mf, &(minint[0])) < 0) ||
             (xtc_int(mf, &(minint[1])) < 0) ||
             (xtc_int(mf, &(minint[2])) < 0) ||
             (xtc_int(mf, &(maxint[0])) < 0) ||
             (xtc_int(mf, &(maxint[1])) < 0) ||
             (xtc_int(mf, &(maxint[2])) < 0))
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return -1;
        }

        sizeint[0] = maxint[0] - minint[0]+1;
        sizeint[1] = maxint[1] - minint[1]+1;
        sizeint[2] = maxint[2] - minint[2]+1;

        /* check if one of the sizes is to big to be multiplied */
        if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
        {
            bitsizeint[0] = xtc_sizeofint(sizeint[0]);
            bitsizeint[1] = xtc_sizeofint(sizeint[1]);
            bitsizeint[2] = xtc_sizeofint(sizeint[2]);
            bitsize       = 0; /* flag the use of large sizes */
        }
        else
        {
            bitsize = xtc_sizeofints(3, sizeint);
        }

        if (xtc_int(mf, &smallidx) < 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return -1;
        }

        maxidx       = MIN(LASTIDX, smallidx + 8);
        minidx       = maxidx - 8; /* often this equal smallidx */
        smaller      = xtc_magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
        smallnum     = xtc_magicints[smallidx] / 2;
        sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx];
        larger       = xtc_magicints[maxidx];

        /* buf[0] holds the length in bytes */

        if (xtc_int(mf, &(buf[0])) < 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return -1;
        }


        if (xtc_data(mf, (char *)&(buf[3]), (int)buf[0]) < 0)
        {
            if (we_should_free)
            {
                free(ip);
                free(buf);
            }
            return -1;
        }



        buf[0] = buf[1] = buf[2] = 0;

        lfp           = fp;
        inv_precision = 1.0 / *precision;
        run           = 0;
        i             = 0;
        lip           = ip;
        while (i < lsize)
        {
            thiscoord = (int *)(lip) + i * 3;

            if (bitsize == 0)
            {
                thiscoord[0] = xtc_receivebits(buf, bitsizeint[0]);
                thiscoord[1] = xtc_receivebits(buf, bitsizeint[1]);
                thiscoord[2] = xtc_receivebits(buf, bitsizeint[2]);
            }
            else
            {
                xtc_receiveints(buf, 3, bitsize, sizeint, thiscoord);
            }

            i++;
            thiscoord[0] += minint[0];
            thiscoord[1] += minint[1];
            thiscoord[2] += minint[2];

            prevcoord[0] = thiscoord[0];
            prevcoord[1] = thiscoord[1];
            prevcoord[2] = thiscoord[2];


            flag       = xtc_receivebits(buf, 1);
            is_smaller = 0;
            if (flag == 1)
            {
                run        = xtc_receivebits(buf, 5);
                is_smaller = run % 3;
                run       -= is_smaller;
                is_smaller--;
            }
            if (run > 0)
            {
                thiscoord += 3;
                for (k = 0; k < run; k += 3)
                {
                    xtc_receiveints(buf, 3, smallidx, sizesmall, thiscoord);
                    i++;
                    thiscoord[0] += prevcoord[0] - smallnum;
                    thiscoord[1] += prevcoord[1] - smallnum;
                    thiscoord[2] += prevcoord[2] - smallnum;
                    if (k == 0)
                    {
                        /* interchange first with second atom for better
                         * compression of water molecules
                         */
                        tmp          = thiscoord[0]; thiscoord[0] = prevcoord[0];
                        prevcoord[0] = tmp;
                        tmp          = thiscoord[1]; thiscoord[1] = prevcoord[1];
                        prevcoord[1] = tmp;
                        tmp          = thiscoord[2]; thiscoord[2] = prevcoord[2];
                        prevcoord[2] = tmp;
                        *lfp++       = prevcoord[0] * inv_precision;
                        *lfp++       = prevcoord[1] * inv_precision;
                        *lfp++       = prevcoord[2] * inv_precision;
                    }
                    else
                    {
                        prevcoord[0] = thiscoord[0];
                        prevcoord[1] = thiscoord[1];
                        prevcoord[2] = thiscoord[2];
                    }
                    *lfp++ = thiscoord[0] * inv_precision;
                    *lfp++ = thiscoord[1] * inv_precision;
                    *lfp++ = thiscoord[2] * inv_precision;
                }
            }
            else
            {
                *lfp++ = thiscoord[0] * inv_precision;
                *lfp++ = thiscoord[1] * inv_precision;
                *lfp++ = thiscoord[2] * inv_precision;
            }
            smallidx += is_smaller;
            if (is_smaller < 0)
            {
                smallnum = smaller;
                if (smallidx > FIRSTIDX)
                {
                    smaller = xtc_magicints[smallidx - 1] /2;
                }
                else
                {
                    smaller = 0;
                }
            }
            else if (is_smaller > 0)
            {
                smaller  = smallnum;
                smallnum = xtc_magicints[smallidx] / 2;
            }
            sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx];
        }
    }
    if (we_should_free)
    {
        free(ip);
        free(buf);
    }
    return mdio_seterror(MDIO_SUCCESS);
}



// other functions to get more info
int xtc_get_water_no(FILE *file_ptr, int *water_index, int water_length, int *wn_index)
{
	fseek(file_ptr, 52, SEEK_CUR);
	int water_sub = 0, *wn_tmp_index;
	wn_tmp_index = (int*)malloc(water_length * sizeof(int));
	


	md_file *mf;
	float *fp, *precision;
	int *size;
	size = (int*)malloc(sizeof(int));
	precision = (float*)malloc(sizeof(float));
	mf = (md_file *)malloc(sizeof(md_file));
	mf->f = file_ptr;
	mf->fmt = MDFMT_XTC;
	mf->mode = XTC_READ;


	int *ip = NULL;
	int oldsize;
	int *buf;

	int minint[3], maxint[3], *lip;
	int smallidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
	int flag, k;
	int small, smaller, i, is_smaller, run;
	float *lfp;
	int tmp, *thiscoord,  prevcoord[3];

	int bufsize, lsize;
	unsigned int bitsize;
	float inv_precision;

    /* avoid uninitialized data compiler warnings */
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;


	if (xtc_int(mf, &lsize) < 0) 
	{
		return -1;
	}

	// if (*size != 0 && lsize != *size) return mdio_seterror(MDIO_BADFORMAT);
	*size = lsize;
	fp = (float*)malloc((lsize * 3) * sizeof(float));

	// *size = lsize;
	size3 = *size * 3;
	if (*size <= 9) {
		for (i = 0; i < *size; i++) {
			if (xtc_float(mf, fp + (3 * i)) < 0) return -1;
			if (xtc_float(mf, fp + (3 * i) + 1) < 0) return -1;
			if (xtc_float(mf, fp + (3 * i) + 2) < 0) return -1;
		}
		return *size;
	}


	xtc_float(mf, precision);
	if (ip == NULL) {
		ip = (int *)malloc(size3 * sizeof(*ip));
		if (ip == NULL) return mdio_seterror(MDIO_BADMALLOC);
		bufsize = (int) (size3 * 1.2);
		buf = (int *)malloc(bufsize * sizeof(*buf));
		if (buf == NULL) return mdio_seterror(MDIO_BADMALLOC);
		oldsize = *size;
	} else if (*size > oldsize) {
		ip = (int *)realloc(ip, size3 * sizeof(*ip));
		if (ip == NULL) return mdio_seterror(MDIO_BADMALLOC);
		bufsize = (int) (size3 * 1.2);
		buf = (int *)realloc(buf, bufsize * sizeof(*buf));
		if (buf == NULL) return mdio_seterror(MDIO_BADMALLOC);
		oldsize = *size;
	}
	buf[0] = buf[1] = buf[2] = 0;

	xtc_int(mf, &(minint[0]));
	xtc_int(mf, &(minint[1]));
	xtc_int(mf, &(minint[2]));

	xtc_int(mf, &(maxint[0]));
	xtc_int(mf, &(maxint[1]));
	xtc_int(mf, &(maxint[2]));
		
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
		bitsizeint[0] = xtc_sizeofint(sizeint[0]);
		bitsizeint[1] = xtc_sizeofint(sizeint[1]);
		bitsizeint[2] = xtc_sizeofint(sizeint[2]);
		bitsize = 0; /* flag the use of large sizes */
	} else {
		bitsize = xtc_sizeofints(3, sizeint);
	}

	xtc_int(mf, &smallidx);
	smaller = xtc_magicints[FIRSTIDX > smallidx - 1 ? FIRSTIDX : smallidx - 1] / 2;
	small = xtc_magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx];

	/* check for zero values that would yield corrupted data */
	if ( !sizesmall[0] || !sizesmall[1] || !sizesmall[2] ) {
		printf("XTC corrupted, sizesmall==0 (case 1)\n");
		return -1;
	}


	/* buf[0] holds the length in bytes */
	if (xtc_int(mf, &(buf[0])) < 0) 
	{
		return -1;
	}

	if (xtc_data(mf, (char *) &buf[3], (int) buf[0]) < 0) 
	{
		return -1;
	}

	buf[0] = buf[1] = buf[2] = 0;

	lfp = fp;
	inv_precision = 1.0f / (*precision);
	run = 0;
	i = 0;
	lip = ip;
	water_sub = 0;
	while (i < lsize) {
		// bing...
		// if(water_sub < water_length && water_index[water_sub] == (i+1))
		// {
		// 	wn_tmp_index[water_sub] = buf[0];
		// 	++water_sub;
		// }


		thiscoord = (int *)(lip) + i * 3;

		if (bitsize == 0) {
			thiscoord[0] = xtc_receivebits(buf, bitsizeint[0]);
			thiscoord[1] = xtc_receivebits(buf, bitsizeint[1]);
			thiscoord[2] = xtc_receivebits(buf, bitsizeint[2]);
		} else {
			xtc_receiveints(buf, 3, bitsize, sizeint, thiscoord);
		}

		i++;
		// bing...
		if(water_sub < water_length && water_index[water_sub] == (i+1))
		{
			wn_tmp_index[water_sub] = buf[0];
			++water_sub;
		}

		thiscoord[0] += minint[0];
		thiscoord[1] += minint[1];
		thiscoord[2] += minint[2];

		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];
 

		flag = xtc_receivebits(buf, 1);
		is_smaller = 0;
		if (flag == 1) {
			run = xtc_receivebits(buf, 5);
			is_smaller = run % 3;
			run -= is_smaller;
			is_smaller--;
		}
		if (run > 0) {
			thiscoord += 3;
			for (k = 0; k < run; k+=3) {
				xtc_receiveints(buf, 3, smallidx, sizesmall, thiscoord);
				i++;
				// bing..
				if(water_sub < water_length && water_index[water_sub] == (i+1))
				{
					wn_tmp_index[water_sub] = buf[0];
					++water_sub;
				}
				thiscoord[0] += prevcoord[0] - small;
				thiscoord[1] += prevcoord[1] - small;
				thiscoord[2] += prevcoord[2] - small;
				if (k == 0) {
					/* interchange first with second atom for better
					 * compression of water molecules
					 */
					tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
					prevcoord[0] = tmp;
					tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
					prevcoord[1] = tmp;
					tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
					prevcoord[2] = tmp;
					*lfp++ = prevcoord[0] * inv_precision;
					*lfp++ = prevcoord[1] * inv_precision;
					*lfp++ = prevcoord[2] * inv_precision;

					if ( !sizesmall[0] || !sizesmall[1] || !sizesmall[2] ) {
						printf("XTC corrupted, sizesmall==0 (case 2)\n");
						return -1;
					}

				} else {
					prevcoord[0] = thiscoord[0];
					prevcoord[1] = thiscoord[1];
					prevcoord[2] = thiscoord[2];
				}
				*lfp++ = thiscoord[0] * inv_precision;
				*lfp++ = thiscoord[1] * inv_precision;
				*lfp++ = thiscoord[2] * inv_precision;
			}
		} else {
			*lfp++ = thiscoord[0] * inv_precision;
			*lfp++ = thiscoord[1] * inv_precision;
			*lfp++ = thiscoord[2] * inv_precision;		
		}
		smallidx += is_smaller;
		if (is_smaller < 0) {
			small = smaller;
			if (smallidx > FIRSTIDX) {
				smaller = xtc_magicints[smallidx - 1] /2;
			} else {
				smaller = 0;
			}
		} else if (is_smaller > 0) {
			smaller = small;
			small = xtc_magicints[smallidx] / 2;
		}
		sizesmall[0] = sizesmall[1] = sizesmall[2] = xtc_magicints[smallidx] ;
	}

	// printf("water_sub: %d\n", water_sub);
	// *wn_index = (int*)malloc((water_length) * sizeof(int));
	memcpy((char*)wn_index, (char*)wn_tmp_index, water_length*sizeof(int));
	// printf("%d  %d\n", wn_index[0], wn_index[1]);
	
	free(wn_tmp_index);
	free(size);
	free(mf);
	free(precision);
	free(fp);
	free(ip);
	free(buf);

	return 1;


}





// xtc_write_uc_frame() - write uncompressed frame info to a file;
int xtc_write_uc_frame(md_file *mf, md_ts *ts) {
	mf->mode = XTC_WRITE;

	int n;
	float f;

	int size = ts->natoms; // explicitly initialized to zero.
	float precision = ts->precision;

	if (!mf || !ts) return mdio_seterror(MDIO_BADPARAMS);
	if (!mf->f) return mdio_seterror(MDIO_BADPARAMS);
	if (mf->fmt != MDFMT_XTC) return mdio_seterror(MDIO_WRONGFORMAT);

	// write magic number
	n = XTC_MAGIC;
	if (xtc_int(mf, &n) < 0) return -1;

	// write number of atoms
	n = ts->natoms;
	if (xtc_int(mf, &n) < 0) return -1;

	// write the simulation step
	n = ts->step;
	if (xtc_int(mf, &n) < 0) return -1;

	// write the time value
	f = ts->time;
	if (xtc_float(mf, &f) < 0) return -1;

	// write the basis vectors of the box
  	if ( (xtc_float(mf, &ts->matri_box[0][0]) < 0) ||
         (xtc_float(mf, &ts->matri_box[0][1]) < 0) ||
         (xtc_float(mf, &ts->matri_box[0][2]) < 0) ||
         (xtc_float(mf, &ts->matri_box[1][0]) < 0) ||
         (xtc_float(mf, &ts->matri_box[1][1]) < 0) ||
         (xtc_float(mf, &ts->matri_box[1][2]) < 0) ||
         (xtc_float(mf, &ts->matri_box[2][0]) < 0) ||
         (xtc_float(mf, &ts->matri_box[2][1]) < 0) ||
         (xtc_float(mf, &ts->matri_box[2][2]) < 0) )
    	return -1;

    /* Now we're left with the job of scaling... */
	for (n = 0; n < ts->natoms * 3; n++)
		ts->pos[n] /= ANGS_PER_NM;

	if (!ts->pos) return mdio_seterror(MDIO_BADPARAMS);
	n = xtc_uc_3dfcoord(mf, ts->pos, &size, &precision);
	if (n < 0) return -1;

	return mdio_seterror(MDIO_SUCCESS);
}



int xtc_uc_3dfcoord(md_file *mf, float *fp, int *size, float *precision)
{
    int     *ip  = NULL;
    int     *buf = NULL;
    int 	 bRead;

    /* preallocate a small buffer and ip on the stack - if we need more
       we can always malloc(). This is faster for small values of size: */
    unsigned     prealloc_size = 3*16;
    int          prealloc_ip[3*16], prealloc_buf[3*20];
    int          we_should_free = 0;

    int          minint[3], maxint[3], mindiff, *lip, diff;
    int          lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int          minidx, maxidx;
    unsigned     sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
    int          flag, k;
    int          smallnum, smaller, larger, i, is_small, is_smaller, run, prevrun;
    float       *lfp, lf;
    int          tmp, *thiscoord,  prevcoord[3];
    unsigned int tmpcoord[30];

    int          bufsize, xdrid, lsize;
    unsigned int bitsize;
    float        inv_precision;
    int          errval = 1;
    int          rc;

    bRead         = (mf->mode == XTC_READ);
    bitsizeint[0] = bitsizeint[1] = bitsizeint[2] = 0;
    prevcoord[0]  = prevcoord[1]  = prevcoord[2]  = 0;

    if (!bRead)
    {
        /* mf is open for writing */
    	// write atoms
        if(xtc_int(mf, size) < 0)
        {
        	return -1;
        }
        size3 = *size * 3;
        /* when the number of coordinates is small, don't try to compress; just
         * write them as floats using xdr_vector
         */
        if (*size <= 9)
        {
            for (i = 0; i < *size; i++) {
				if (xtc_float(mf, fp + (3 * i)) < 0) return -1;
				if (xtc_float(mf, fp + (3 * i) + 1) < 0) return -1;
				if (xtc_float(mf, fp + (3 * i) + 2) < 0) return -1;
			}
			return mdio_seterror(MDIO_SUCCESS);
        }

        int check_num = 7750;
        if(fwrite((char*)&check_num, 4, 1, mf->f) != 1){
        	printf("%s %d: check_num write failed!\n", __FUNCTION__, __LINE__);
        }


        int uc_len = size3 * sizeof(float);
        if (xtc_uc_data(mf, (char*)fp, uc_len) < 0) return -1;
        
        return mdio_seterror(MDIO_SUCCESS);

    }
    else
    {
    	return mdio_seterror(MDIO_MODEERROR);
    }

    return mdio_seterror(MDIO_SUCCESS);
}


int xtc_uc_data(md_file *mf, char *buf, int len)
{
	if (!mf || !buf || len < 1) return mdio_seterror(MDIO_BADPARAMS);
	size_t slen = (size_t)len;

	if (fwrite(buf, 1, slen, mf->f) != slen) {
		if (feof(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_EOF);
		}
		else if (ferror(mf->f)) {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_IOERROR);
		}
		else {
			XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
			return mdio_seterror(MDIO_UNKNOWNERROR);
		}
	}

	if (len % 4) {
		// round bytes count to full xtc/xdr units
		int rndup = BYTES_PER_XTC_UNIT - (len % 4);
		if (fwrite(xtc_zero, rndup, 1, mf->f) != 1) {
			if (feof(mf->f)) {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_EOF);
			}
			else if (ferror(mf->f)) {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_IOERROR);
			}
			else {
				XTC_DEBUG(__FILE__, __FUNCTION__, __LINE__);
				return mdio_seterror(MDIO_UNKNOWNERROR);
			}
		}
	}
	return len;
}