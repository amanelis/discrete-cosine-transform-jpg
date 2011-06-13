/*
Written by: Alex Manelis
Subject: For Professor Ortiz's
CS 4953 Steganography Class
Summer 2009
*/

/* INCLUDE DECLARATIONS */
#include <stdio.h>
#include <math.h>

/* VARIABLE DECLARATIONS */
static const double PI = 3.141593;
static const int COL = 8;
static const int ROW = 8;

double dctblock[COL][ROW];
double idctblock[COL][ROW];
double result[8][8][8][8];

/****** NEW ASSIGNMENT******/
unsigned char uc_block[COL][ROW];
double dc_block[COL][ROW];
unsigned char idct_block[COL][ROW];

void init_rand();
void init_zero();
void init_twff();
void print_uc_block();
void print_dc_block();
void print_idct_block();
void getDiff(unsigned char buc_array[][8], unsigned char idc_array[][8]);
unsigned char dtuc(double n);
/***************************/

/* FUNCTION DECLARATIONS */
void printblock(double values[][8]);
void dct(double values[][8]);
void idct(double values[][8]);
void basis();
double calCos(double a, double b);
double c(int);
void diffMtx(double oblock[][8], double iblock[][8]);

/* MAIN */
int main(){    
    //Pseudo random valued array
	double ablock[COL][ROW] = {
										{140,144,147,140,140,155,179,175},
										{144,152,140,147,140,148,167,179},
										{152,155,136,167,163,162,152,172},
										{168,145,156,160,152,155,136,160},
										{162,148,156,148,140,136,147,162},
										{147,167,140,155,155,140,136,162},
										{136,156,123,167,162,144,140,147},
										{148,155,136,155,152,147,147,136}
    };
    
    //Pseudo random valued array
	double bblock[COL][ROW] = {
										{16,11,10,16,24,40,51,61},
										{12,12,14,19,26,58,60,55},
										{14,13,16,24,40,57,69,56},
										{14,17,22,29,51,87,80,62},
										{18,22,37,56,68,109,103,77},
										{24,35,55,64,81,104,113,92},
										{49,64,78,87,103,121,120,101},
										{72,92,95,98,112,100,103,90}
	};
    
    //All zero array
    double cblock[COL][ROW] = {
    
                                        {0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0}
    };
    
    //All 255 array
    double dblock[COL][ROW] = {

                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255}
    };

    //Array with +255 values
    double eblock[COL][ROW] = {

                                        {289,345,269,453,299,645,334,356},
                                        {889,945,669,853,399,945,534,956},
                                        {489,845,769,653,599,445,534,956},
                                        {689,745,369,353,999,545,734,856},
                                        {389,545,469,953,299,345,834,356},
                                        {389,945,869,853,599,545,934,756},
                                        {289,545,969,653,499,845,334,756},
                                        {189,445,369,553,799,745,434,756}
    };

    //Unsigned char array
    unsigned char carray[COL][ROW] = {
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255},
                                        {255,255,255,255,255,255,255,255}
    };   
   
    //Initializes the starting array 
    init_rand();

    //Calculates the basis 4x4 array values
    basis();

    //Prints out the starting block that will be used
    printf("Starting Block:\n");
    print_uc_block();

    //Calculate DCT and store value into the GLOBAL dc_block[][]
    for(int u=0;u<8;u++) {
     		for(int v=0;v<8;v++) {
                    double temp = 0.0;
                    for(int x=0;x<8;x++) {
                            for(int y=0;y<8;y++) {
                                    temp += (uc_block[x][y] - 128) * result[u][v][x][y];
                            }
                    }
                    dc_block[u][v] = temp * (c(u)) * (c(v)) * .25;
            }
    }
    
    //Print the DCT BLOCK
    printf("DCT Block:\n");
    print_dc_block();

    //Caculate the IDCT and write values back into idct_block[][]
    for(int x=0;x<8;x++) {
            for(int y=0;y<8;y++) {
                    double temp = 0.0;
                    for(int u=0;u<8;u++) {
                            for(int v=0;v<8;v++) {
                                temp += (c(u)) * (c(v)) * dc_block[u][v] * cos(PI * (2 * x + 1) * u / 16) * cos(PI * (2 * y + 1) * v / 16);
                            }
                    }
                    //This round function is what makes the 
                    //Differnece matrix output perfectly matching
                    //values, if not used there will be at most (1)
                    //value of difference.
                    idct_block[x][y] = round((.25 * temp) + 128);
                    
                    //This assignment's difference matrix will output
                    //at most (1) at about 50% of the matrix's values
                    //idct_block[x][y] = ((temp * .25) + 128);
            }
    }

    //Print out original from IDCT
    printf("IDCT Block:\n");
    print_idct_block();

    //Do conversions and compute the difference matrix
    printf("Difference Matrix:\n");
    getDiff(uc_block, idct_block);
    
    /*
    //basis() calculates base values that will be used in idct function
    //to caluclate the inverse of our dctblock array
    basis();

    //Standard array print function for output   
    printf("Starting block:\n"); 
    printblock(ablock);

    //dct(block) calculates dct and then outputs to global array dctblock
    dct(ablock);

    //Print output of array
    printf("DCT block:\n");
    printblock(dctblock);
 
    //takes the dctblock from the previous function and calculates
    //the inverse and uses the values stored in result[][][][]
    //then output result for idct into the global idctblock
    idct(dctblock);

    //Print output, should be the same as initial block on first print
    printf("IDCT block:\n");
    printblock(idctblock);

    printf("Difference of Starting and IDCT\n");
    diffMtx(ablock, idctblock);
    */
	return 0;
}
//Does conversion for double -> unsigned char
//NOT ACTUALLY USED AFTER I REALIZED NO CONVERSION NEEDED
unsigned char dtuc(double n) {
        int i = (int)n;
        return (i + '0' + ((n-i) > 0.4999));
}


//Finds the difference in the two arrays, NOT NEEDED -> also does conversion double->unsigned char
void getDiff(unsigned char buc_array[][8], unsigned char idc_array[][8]){
    unsigned char diffOut[COL][ROW];
    unsigned char temp_buc;  //Holds value of usigned char array
    unsigned char temp_idc;  //Holds value of dct->unsigned char array
    unsigned char output;
    
    for(int a=0; a<COL; a++){
        for(int b=0; b<ROW; b++){
            temp_buc, temp_idc, output = 0;
            temp_buc = buc_array[a][b];
            temp_idc = idc_array[a][b];
            output = temp_buc - temp_idc;
            printf("%d\t", output);
        }
        printf("\n");
    }
}

//Printing the unsigned char block
void print_uc_block(){
    int i, j;
    for(i=0; i<COL; i++){
        for(j=0; j<ROW; j++){
            printf("%d\t", uc_block[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

//Printing the int dct block
void print_dc_block(){
    int i, j;
    for(i=0; i<COL; i++){
        for(j=0; j<ROW; j++){
            printf("%6.1f\t", dc_block[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void print_idct_block(){
    int i, j;
    for(i=0; i<COL; i++){
        for(j=0; j<ROW; j++){
            printf("%d\t", idct_block[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

//Standard array printing, with a little formatting to distinguish col, row
void printblock(double values[][8]) {
	int a, b;	
	for(a=0; a<COL; a++){
		for(b=0; b<ROW; b++){
            //printf("%6.1f\t", values[a][b]);
            printf("%6.1f  ", values[a][b]);
            if(b == 7){
                printf("\n");
            }
		}
	}
	printf("\n");
}

//Input is array in main
//Writes out to GLOBAL array: dctblock[8][8]
void dct(double values[][8]) {
	int u, v, x, y;
	for(u=0;u<8;u++) {
		for(v=0;v<8;v++) {
		    double temp = 0.0;
			for(x=0;x<8;x++) {
				for(y=0;y<8;y++) {
					//temp += values[x][y] * cos(PI *(2 * x + 1) * u /16 ) * cos(PI * (2 * y + 1) * v / 16);
					temp += values[x][y] * result[u][v][x][y];
				}
			}
            //store values in GLOBAL array dctblock
			dctblock[u][v] = temp * (c(u)) * (c(v)) * 1/4;
		}
	}
	printf("\n");
}

//Input is the GLOBAL array: dctblock[8][8]
//Writes out to GLOBAL array: idctblock[8][8]
void idct(double values[][8]) {
	int u, v, x, y;
	for(x=0;x<8;x++) {
		for(y=0;y<8;y++) {
			double temp = 0.0;
			for(u=0;u<8;u++) {
				for(v=0;v<8;v++) {
                    //using already pre generated values in the result array
				    temp += (c(u)) * (c(v)) * values[u][v] * result[u][v][x][y];
                }
			}
            //store values in GLOBAL array idctblock
			idctblock[x][y] = temp * 1/4;
		}
	}
    printf("\n");
}

//basis() caluclates the base values on the array indicies
//used in the inverse valuculation
void basis() {
    int u, v, x, y;
    double temp = 0;
    for(u=0; u<8; u++){
        for(v=0; v<8; v++){
            for(x=0; x<8; x++){
                for(y=0; y<8; y++){
                    //sends values into calCos to get half of the equation completed
                    //then multiplies result and puts it into result[][][][]
                    result[u][v][x][y] = calCos(x, u) * calCos(y, v);
                    /*
                    printf("%6.1f\t", result[u][v][x][y]);
                    if(y == 7)
                        printf("\n");
                    */
                }
            }
        }
    }
}

//calCos, function to calculate the cos of one side of the cos equation
double calCos(double a, double b){
    double temp = 0;
    temp = (double) cos(PI * (2 * a + 1) * b / 16);
    return temp;
}

//c() handles the 0 or 1 case of the equation
double c(int number) {
	if(number == 0) {
		return 1/sqrt(2);
	} else {
		return 1;
	}
}

//matrix difference
void diffMtx(double oblock[][8], double iblock[][8]){
    int a, b, c, d;
    double orig, invdct;
    for(a=0; a<8; a++){
        for(b=0; b<8; b++){
            orig, invdct = 0.0;
            double o, i = 0.0;
            o = oblock[a][b];
            i = iblock[a][b];
            
            if(o == i){
                printf("0  ");
            } else {
                printf("1  ");
            }

            if(b == 7)
                printf("\n");
            
        }
    }
    printf("\n");
}

//Initializes the starting block with pseudo random values
void init_rand(){
    //UHHHHHHHHHHHH why wont this assignment work.......?!?!
    /*
    uc_block[COL][ROW] = {
										{140,144,147,140,140,155,179,175},
										{144,152,140,147,140,148,167,179},
										{152,155,136,167,163,162,152,172},
										{168,145,156,160,152,155,136,160},
										{162,148,156,148,140,136,147,162},
										{147,167,140,155,155,140,136,162},
										{136,156,123,167,162,144,140,147},
										{148,155,136,155,152,147,147,136}
    };
    */
	

    uc_block[0][0] = 140;
    uc_block[1][0] = 144;
    uc_block[2][0] = 147;
    uc_block[3][0] = 140;
    uc_block[4][0] = 140;
    uc_block[5][0] = 155;
    uc_block[6][0] = 179;
    uc_block[7][0] = 175;
    uc_block[0][1] = 144;
    uc_block[1][1] = 152;
    uc_block[2][1] = 140;
    uc_block[3][1] = 147;
    uc_block[4][1] = 140;
    uc_block[5][1] = 148;
    uc_block[6][1] = 167;
    uc_block[7][1] = 179;
    uc_block[0][2] = 152;
    uc_block[1][2] = 155;
    uc_block[2][2] = 136;
    uc_block[3][2] = 167;
    uc_block[4][2] = 163;
    uc_block[5][2] = 162;
    uc_block[6][2] = 152;
    uc_block[7][2] = 172;
    uc_block[0][3] = 168;
    uc_block[1][3] = 145;
    uc_block[2][3] = 156;
    uc_block[3][3] = 160;
    uc_block[4][3] = 152;
    uc_block[5][3] = 155;
    uc_block[6][3] = 136;
    uc_block[7][3] = 160;
    uc_block[0][4] = 162;
    uc_block[1][4] = 148;
    uc_block[2][4] = 156;
    uc_block[3][4] = 148;
    uc_block[4][4] = 140;
    uc_block[5][4] = 136;
    uc_block[6][4] = 147;
    uc_block[7][4] = 162;
    uc_block[0][5] = 147;
    uc_block[1][5] = 167;
    uc_block[2][5] = 140;
    uc_block[3][5] = 155;
    uc_block[4][5] = 155;
    uc_block[5][5] = 140;
    uc_block[6][5] = 136;
    uc_block[7][5] = 162;
    uc_block[0][6] = 136;
    uc_block[1][6] = 156;
    uc_block[2][6] = 123;
    uc_block[3][6] = 167;
    uc_block[4][6] = 162;
    uc_block[5][6] = 144;
    uc_block[6][6] = 140;
    uc_block[7][6] = 147;
    uc_block[0][7] = 148;
    uc_block[1][7] = 155;
    uc_block[2][7] = 136;
    uc_block[3][7] = 155;
    uc_block[4][7] = 152;
    uc_block[5][7] = 147;
    uc_block[6][7] = 147;
    uc_block[7][7] = 136;
}

//Initializes the starting block with all zero's
void init_zero(){
    for(int a=0; a<COL; a++){
		for(int b=0; b<ROW; b++){
			uc_block[a][b] = 0;
		}
	}
}

//Initializes the starting block with all 255's
void init_twff(){
    for(int a=0; a<COL; a++){
		for(int b=0; b<ROW; b++){
			uc_block[a][b] = 255;
		}
	}
}

