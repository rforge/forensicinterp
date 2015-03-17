#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>

/*
Connectiviy measure of a pixel in a matrix

It Computes the connectiviy measure of a pixel in a matrix of 0's and 1's.
The connectiviy measure is taken from Yokoi et al.

Parameters:
- 'matrix' A numeric matrix, representing the black and white image (a pixel of value 0 is from the background).
- 'i' A numeric value, representing the index of the pixel to analyse.
- 'w' A numeric value, the image's width.
- 'h' A numeric value, the image's height.

Return:
The connectivity measure C of the pixel based on Yokoi et al.

Reference:
S. Yokoi, J. Toriwaki, T. Fukurama  (1973), "Topological properties in digitized binary pictures." Systems Computers Controls, 4, pp. 32-39.

Author:
Alexandre Thiery
*/
int GetConnectivity(int matrix[], int i , int w, int h)
{
    return (matrix[i + h] - (matrix[i + h] * matrix[i - 1 + h] * matrix[i - 1])) +
        (matrix[i - 1] - (matrix[i - 1] * matrix[i - 1 - h] * matrix[i - h])) +
        (matrix[i - h] - (matrix[i - h] * matrix[i + 1 - h] * matrix[i + 1])) +
        (matrix[i + 1] - (matrix[i + 1] * matrix[i + 1 + h] * matrix[i + h]));
}


/*
Extract the skeleton of a black and white image.
 
The skeletonisation algorithm implemented is based on Stentiford (1983). It is a recursive thinning process that allows to extract a curve with a thickness of 1 pixel.
In the original version of this algorithm, endpoints pixels 1 could not be removed (endpoint pixel are defined as a black element that are 8–connected to only one black element). The algorithm used in this package enables the removal of such pixels, in order to extract closed loops only.

@param matrix A numeric matrix, representing the black and white image (a pixel of value 0 is from the background).
@return A matrix with values of 1 for the contour and values of 0 for the background. 

Reference:
F.W.M. Stentiford, R.G. Mortimer (1983), "Some new heuristics for thinning binary handprinted characters for OCR." IEEE Transactions on Systems, Man and Cybernetics, 13, pp. 81-84.

Author:
Alexandre Thiery
*/

void ImageThin(int *matrix, int *nrow, int *ncol) 
{
	int nr = *nrow;
	int nc = *ncol;

	int change = 1;
	//int *destination = malloc(sizeof(int)*nc*nr);
	int destination[nr*nc];
	//int matrix[nr*nc];
	//for(int i = 0; i < (nr*nc); ++i)
	//	matrix[i] = source[i];
	for(int i = 0; i < (nr*nc); ++i)
		destination[i] = matrix[i];

	int index;
	int count = 0;

	//if(verbose == 1) Rprintf("param: nr = %d nc=%d\n", nr, nc);
	while(change == 1)
	{
		//if(verbose == 1) Rprintf("begin\n");
		change = 0;

		index = 0;
		for(int i = 0; i < nr; ++i) {
		for(int j = 0; j < nc; ++j) {
			if(destination[index] == 1 & (matrix[index + 1] == 1 & matrix[index - 1] == 0))
			{
				if(GetConnectivity(destination, index, nc, nr) == 1)
				{
					destination[index] = 0;
					change = 1;	
					//Rprintf("%d %d (%d  %d %d)\n", destination[index], matrix[index] , i,j,index);
				}
			}
			index++;
		}
		}


		index = 0;
		for(int i = 0; i < nr; ++i) {
		for(int j = 0; j < nc; ++j, ++index) {
			if(destination[index] == 1 & (matrix[index + nr] == 0 & matrix[index - nr] == 1))
			{
				if(GetConnectivity(destination, index, nc, nr) == 1)
				{
					destination[index] = 0;
					change = 1;	
				}
			}
		}
		}

		index = 0;
		for(int i = 0; i < nr; ++i) {
		for(int j = 0; j < nc; ++j, ++index) {
			if(destination[index] == 1 & (matrix[index + nr] == 1 & matrix[index - nr] == 0))
			{
				if(GetConnectivity(destination, index, nc, nr) == 1)
				{
					destination[index] = 0;
					change = 1;	
	
				}
			}
		}
		}


		index = 0;
		for(int i = 0; i < nr; ++i) {
		for(int j = 0; j < nc; ++j, ++index) {
			if(destination[index] == 1 & (matrix[index + 1] == 0 & matrix[index - 1] == 1))
			{
				if(GetConnectivity(destination, index, nc, nr) == 1)
				{
					destination[index] = 0;
					change = 1;	
				}
			}
		}
		}

	
		for(int i = 0; i < (nr*nc); ++i)
			matrix[i] = destination[i];

		count++;

	//if(verbose == 1) Rprintf("end %d\n",count);
	//if(count > 10) change = 0;

	}

	/*for(int i = 0; i < (nr*nc); ++i)
	{
		if(destination[i] == 1)
			source[i] = 1;
		else source[i] = 0;		
	}*/
	//matrix[i] = destination[i];
}
