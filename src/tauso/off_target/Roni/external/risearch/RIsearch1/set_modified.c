
#include "set_modified.h"

char buffer[32000];
extern short modified_dsm[6][6][6][6];
extern short modified_dsm_extend[6][6][6][6];
short* p_modified_dsm = &modified_dsm[0][0][0][0]; 
short* p_modified_dsm_extend = &modified_dsm_extend[0][0][0][0];

  
void set_modified_dsm(short* p_dsm,const char* csv_filename){
	    FILE *file;
	    short number;
	    file = fopen( csv_filename , "r");
	    if (file == NULL) {
	        printf("Could not open the file.\n");
	    }
	    int i=0;
	    while (fscanf(file, "%hd,", &number) != EOF) {
	        *(p_dsm+i)=number;
	        i++;
	    }
	    fclose(file);
}

void set_modified() {
		set_modified_dsm(p_modified_dsm,"./modified_dsm.csv");
		set_modified_dsm(p_modified_dsm_extend,"./modified_dsm_extend.csv");
}

void test(){
	printf("running test:\n");
	extern short dsm_su95_rev_woGU_pos[6][6][6][6];
	short* p_dsm_su95_rev_woGU_pos =&dsm_su95_rev_woGU_pos[0][0][0][0];
	set_modified();
	int counter = 0 ;
	for(int i=0 ; i<1296 ;i++){
		int x = *(p_modified_dsm +i);
		int y =  *(p_dsm_su95_rev_woGU_pos+i);
		if( x != y ) printf("%d != %d index: %d\n",x,y,i);
		else counter++;	
	}
	if(counter == 1296) printf("PASS\n");
	else printf("FAILED\n");	
}
