#include "plfs_xtc.h"
#include <stdlib.h>
#include <string.h>

// for test;
void print_mdh(md_header *mdh);
void print_mdts(md_ts *mt);


int main()
{
	md_file *mf, *mf_write, *mf_read, *uc_mf;
	md_header mdh;
	md_ts mt;

	// test read frame 
	printf("----> test read frmae <----\n");
	mf = xtc_open("~/Data-for-VMD/water-channel/step7_extend-400ns-fit-390ns-allwater.xtc", MDFMT_XTC, XTC_READ);
	uc_mf = xtc_open("/tmp/temp.xtc", MDFMT_XTC, XTC_WRITE);

	int count = 0;
	while(1)
	{
		if(feof(mf->f)){
			printf("EOF  %d\n", __LINE__);
			break;
		}
		if(xtc_read_frame(mf, &mt))
		{
			printf("%s\n", mdio_errmsg(mdio_errno()));
			if(feof(mf->f)){
				printf("EOF  %d\n", __LINE__);
			}
			break;
		}
		++count;
		printf("%d\n", count);

		// write uncompressed info;
		if(xtc_write_uc_frame(uc_mf, &mt))
		{
			printf("%d: %s\n", __LINE__, mdio_errmsg(mdio_errno()));
			break;
		}

		free(mt.pos);
	}
	
	if(xtc_close(mf) && xtc_close(uc_mf))
	{
		printf("%s\n", mdio_errmsg(mdio_errno()));
		exit(1);
	}	
	
	return 0;
}



void print_mdh(md_header *mdh)
{
	printf("----> md_header: %d  %f\n", mdh->natoms, mdh->timeval);
}

void print_mdts(md_ts *mt)
{
	printf("----> natoms: %d  step: %d  time: %f  precision: %f\n", mt->natoms, mt->step, mt->time, mt->precision);
	printf("----> md_box:\nA: %f  B: %f  C: %f\nalpha: %f  beta: %f  gamma: %f\n",
		mt->box->A, mt->box->B, mt->box->C, 
		mt->box->alpha, mt->box->beta, mt->box->gamma );
	printf("----> pos:\n");
	// for(int i = 0; i < 20; ++i)
	// {
	// 	printf("index: %d ----> %f %f %f\n", i+1, mt->pos[i*3], mt->pos[i*3+1], mt->pos[i*3+2]);
	// }
}
