// mexVolWrite.cc
// DESCRIPTION: writing of binary volume files
// AUTHOR: Daniel Haenschke, ISS, KIT Karlsruhe, Germany
// 29/06/2009: release version
//
// to do:
//  - input type has to be double, should be arbitrary and converted
//    to requested output type
//

#include "mex.h"
#include <string>

struct DataTypeTableStruct {
    std::string key;
    int value;
    short dataLen;
};

enum DataTypes {
    UINT8, FLOAT32
};

static const struct DataTypeTableStruct DataTypeTable [] =
{
    { "uint8",	UINT8,  1 },
    { "float32",FLOAT32,  4 },
    { "endOfStruct", -1, -1 }
};

int lookup_DataTypeTable (const struct DataTypeTableStruct *tbl, std::string& search_str)
{
    int k = -1;
    while (!(tbl[++k].key.compare("endOfStruct")==0)) {
	if (tbl[k].key.compare(search_str)==0)
	    return k;
    }
    return -1; // no match found
}

inline mwSize getPos3D(mwSize xPos, mwSize yPos, mwSize zPos, const mwSize *volDims)
{
 return(zPos*volDims[1]*volDims[0]+yPos*volDims[0]+xPos);
}

void printUsage()
{
 mexPrintf("\nSyntax:\nout = mexVolWrite(filename,vol,datatype)\n\n");
 return;
}

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
 (void) plhs;   // not used
 if (nrhs != 3) {
	printUsage();
	mexErrMsgTxt("3 input arguments required:\n-filename (string)\n-size (3dim vector, unsigned integer)\n-data type (string: 'float32' or 'uint8')");
 }
 if (nlhs > 0) {
	printUsage();
	mexErrMsgTxt("Too many output arguments.\nNo output is given back.");
 }

  if (!(mxIsChar(prhs[2]))){
	printUsage();
        mexErrMsgTxt("Third input argument must be of type string.");
    }
  char *strIn2;
  strIn2 = mxArrayToString(prhs[2]);
  std::string datatype(strIn2);
  mxFree(strIn2);

  if (!(mxIsChar(prhs[0]))){
	printUsage();
        mexErrMsgTxt("First input argument must be of type string.");
    }
  char *strIn;
  strIn = mxArrayToString(prhs[0]);
  std::string filename(strIn);
  mxFree(strIn);

  mwSize inND = mxGetNumberOfDimensions(prhs[1]);
  if ( inND!=3 || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
	printUsage();
	mexErrMsgTxt("Input has to be a real 3dim array of type double.");
  }
  const mwSize *inDim = mxGetDimensions(prhs[1]);
 
  double *prFrom;
  prFrom=(double *)mxGetPr(prhs[1]);
  mwSize x,y,z;

  int k = lookup_DataTypeTable(DataTypeTable,datatype);
  if (k<0) {
	printUsage();
	mexErrMsgTxt("\nerror!\nthe third argument must be a one of the following strings:\n - uint8\n - float32");
  }
  short dataLen = DataTypeTable[k].dataLen;
  int dataType = DataTypeTable[k].value;

  FILE * inp;
  inp = fopen(filename.c_str(), "wb");
  if (inp==NULL) {
	char errmsg[256];
	printUsage();
	sprintf(errmsg,"\ncannot create output file \"%s\"", filename.c_str());
	mexErrMsgTxt(errmsg);
  }

  int dataSize = inDim[0]*dataLen;

  char *buffer = new char [dataSize];

  if (!buffer) {
    	printUsage();
	mexErrMsgTxt("\nnot enough free memory for buffered writing of given volume");
  }

  char *b = buffer;

  for (z = 0; z < inDim[2]; z++) {
    for (y = 0; y < inDim[1]; y++) {
	for (x = 0; x < inDim[0]; x++) {
	    switch (dataType) {
		case UINT8:	*(unsigned char*)b=(unsigned char)prFrom[getPos3D(x,y,z,inDim)]; b+=dataLen; break;
		case FLOAT32:	*(float*)b=(float)prFrom[getPos3D(x,y,z,inDim)]; b+=dataLen; break;
		default: break;
	    }
	}
	b = buffer;
	if (dataSize != fwrite(buffer, 1, dataSize, inp)) {
		delete [] buffer;
		printUsage();
		mexErrMsgTxt("\nthe volume file cannot be written completely");
	}
	
    }
  }
  fclose(inp);
  delete [] buffer;

  return;
}
