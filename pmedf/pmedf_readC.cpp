 /******************** pmedf_readC.cpp *********************
 *
 * Copyright (c) 2006-2008 Petr Mikulik
 *	Department of Condensed Matter Physics,
 *	Masaryk University, Brno
 *	mikulik@physics.muni.cz
 *
 * Version: 11.5.2008
 *
 * Fast C++ reader of edf files for Octave.
 * It can read .edf, .ehf, .edf.gz, .ehf.gz files.
 *
 * History:
 *	- 11.5.2008: Fixed memory leaks. Long datatypes have size of 4 B.
 *	- 5.7.2006: Fixed s.data() => s.c_str() for std::string.
 *	- 28.6.2006: First version.
 *
 *
 * License: GNU GPL
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
*/


#include <oct.h>
#include <string>
#include <zlib.h>

// debugging?
//#define DEBUG 1
#ifdef DEBUG
#  define OK octave_stdout <<"OK at " <<__LINE__ <<"\n";
#endif


/***************************************************************************
 *
 * Helper functions
 *
 ***************************************************************************/

/* Find value_ptr as pointer to the parameter of the given key in the header.
 * Returns NULL on success.
 */
static char*
edf_findInHeader ( const char* header, const char* key )
{
    char *value_ptr = strstr( header, key );
    if (!value_ptr) return NULL;
    /* an edf line is "key     = value ;" */
    value_ptr = 1 + strchr( value_ptr + strlen(key), '=' );
    while (isspace(*value_ptr)) value_ptr++;
    return value_ptr;
}

/* Routine byteswap: needed for big <-> little endian conversion
 */
static inline void
byteswap ( char* data, int datalen )
{
    if (datalen < 2) return;
    char tmp;
    char *dest = data + datalen - 1;
    while (dest > data) {
	tmp = *dest;
	*dest-- = *data;
	*data++ = tmp;
    }
}


/***************************************************************************
 * Data types and constants
 ***************************************************************************/

#define LowByteFirst  0
#define HighByteFirst 1

// Synonyms:
#define LittleEndian LowByteFirst
#define BigEndian    HighByteFirst

/*
 * Determine the endian type of the machine. This is a runtime detection.
 * Supported are Little and Big endian machines, not the Middle ones (e.g. PDP).
 */

static int int_1 = 1;
#define    THIS_ENDIAN_BYTE_ORDER (((char*)&int_1)[0] ? LowByteFirst : HighByteFirst)
// #define THIS_ENDIAN_BYTE_ORDER (((char*)&int_1)[0] ? LittleEndian : BigEndian)



/* List of EDF supported datatypes
 */
enum DataType {
	U_CHAR_DATATYPE = 0,	CHAR_DATATYPE,	//  8 bits = 1 B
	U_SHORT_DATATYPE, 	SHORT_DATATYPE,	// 16 bits = 2 B
	U_INT_DATATYPE,		INT_DATATYPE,	// 24 bits = 3 B
	U_L_INT_DATATYPE,	L_INT_DATATYPE,	// 32 bits = 4 B
	FLOAT_DATATYPE,		DOUBLE_DATATYPE,// 4 B, 8 B
	UNKNOWN_DATATYPE = -1
};

/* Note - compatibility:
	Unsigned8 = 1,	Signed8,	Unsigned16,	Signed16,
	Unsigned32,	Signed32,	Unsigned64,	Signed64,
	FloatIEEE32,	DoubleIEEE64
*/



/***************************************************************************
 * Tables
 ***************************************************************************/

// table key-value structure
struct table {
    const char *key;
    int value;
};

struct table3 {
    const char *key;
    int value;
    short sajzof;
};
 
/* Returns index of the table tbl whose key matches the beginning of the
 * search string search_str.
 * It returns index into the table or -1 if there is no match.
 */
static int
lookup_table_nth ( const struct table *tbl, const char *search_str )
{
    int k = -1;
    while (tbl[++k].key)
	if (tbl[k].key && !strncmp(search_str, tbl[k].key, strlen(tbl[k].key)))
	return k;
    return -1; // not found
}

static int
lookup_table3_nth ( const struct table3 *tbl, const char *search_str )
{
    int k = -1;
    while (tbl[++k].key)
	if (tbl[k].key && !strncmp(search_str, tbl[k].key, strlen(tbl[k].key)))
	    return k;
    return -1; // not found
}


/***************************************************************************
 * Options for tables
 ***************************************************************************/

static const struct table3 edf_datatype_table[] =
{
    { "UnsignedByte",	U_CHAR_DATATYPE,  1 },
    { "SignedByte",	CHAR_DATATYPE,    1 },
    { "UnsignedShort",	U_SHORT_DATATYPE, 2 },
    { "SignedShort",	SHORT_DATATYPE,   2 },
    { "UnsignedInteger",U_INT_DATATYPE,   4 },
    { "SignedInteger",	INT_DATATYPE,	  4 },
    { "UnsignedLong",	U_L_INT_DATATYPE, 4 },
    { "SignedLong",	L_INT_DATATYPE,   4 },
    { "FloatValue",	FLOAT_DATATYPE,   4 },
    { "DoubleValue",	DOUBLE_DATATYPE,  8 },
    { "Float",		FLOAT_DATATYPE,   4 }, // Float and FloatValue are synonyms
    { "Double",		DOUBLE_DATATYPE,  8 }, // Double and DoubleValue are synonyms
    { NULL, -1, -1 }
};


static const struct table edf_byteorder_table[] =
{
    { "LowByteFirst",	LittleEndian }, /* little endian */
    { "HighByteFirst",	BigEndian },    /* big endian */
    { NULL, -1 }
};

/* Orientation of axes of the raster, as the binary matrix is saved in 
 * the file. (Determines the scanning direction, or the "fastest" index
 * of the matrix in the data file.)
 */
enum EdfRasterAxes {
    RASTER_AXES_XrightYdown,	// matricial format: rows, columns
    RASTER_AXES_XrightYup	// cartesian coordinate system
    // other 6 combinations not available (not needed until now)
};

static const struct table rasteraxes_table[] =
{
    { "XrightYdown", RASTER_AXES_XrightYdown },
    { "XrightYup",   RASTER_AXES_XrightYup },
    { NULL, -1 }
};



/*
 *
 * Implementation of the edf reader octave function.
 *
 */

DEFUN_DLD (pmedf_readC, args, nargout,
"\
Read an edf file\n\
\n\
Usage:\n\
  [header {, image}] = pmedf_readC(edf_filename)\n\
\n\
Reading ESRF header files: .ehf, .edf, .edf.gz, .ehf.gz files.\n\
\n\
Usage:\n\
        [header, data] = pmedf_readC('hello.edf');\n\
        header = pmedf_readC('hello.edf');\n\
\n\
The first syntax reads both the header and the data; the second one reads\n\
only the header, which can be useful for parsing header information.\n\
\n\
Author: Petr Mikulik\n\
Version: 11. 5. 2008\n\
\n\
Example: [h,a] = pmedf_readC('demo.edf');\n\
         imagesc(a');\n\
")
{
    octave_value_list retval = octave_value(0);

    int nargin = args.length();

    int k;
    char *header = NULL;
    int header_size = 0;
    char *p;
    gzFile inp;

    if (nargin != 1) {
	error("Usage:\n    [header {, image}] = pmedf_readC(edf_filename)");
	// print_usage("pmedf_readC");
	return retval;
    }
    if (!(args(0).is_string())) {
	error("The argument must be a string");
	return retval;
    }

    std::string filename = args(0).string_value();
    inp = gzopen(filename.c_str(), "rb");
    if (!inp) {
	error("Cannot open input file \"%s\"", filename.c_str());
	return retval;
    }

    // read header: it is a multiple of 512 B ending by "}\n"
    while (header_size == 0 || strncmp(&header[header_size-2],"}\n",2)) {
	int header_size_prev = header_size;
	if (header_size > 12*512) { /* protection against infinite loop */
	    error("Damaged EDF header of %s: not multiple of 512 B.\n", filename.c_str());
	    free(header);
	    return retval;
	}
	header_size += 512;
	if (!header)
	    header = (char*)malloc(header_size+1);
	else
	    header = (char*)realloc(header, header_size+1);
	header[header_size_prev] = 0; /* protection against empty file */
	// fread(header+header_size_prev, 512, 1, fp);
	gzread(inp, header+header_size_prev, 512);
	header[header_size] = 0; /* end of string: protection against strstr later on */
    }

   // parse the header 
   int dim1 = -1, dim2 = -1, datatype = -1, datalen = -1;
   char *otherfile_name = 0; // this file, or another file with the data (EDF vs EHF formats)
   int otherfile_skip = 0;

    if ((p = edf_findInHeader(header, "EDF_BinaryFileName"))) {
	int plen = strcspn(p, " ;\n");
	otherfile_name = (char*)realloc(otherfile_name, plen+1);
	strncpy(otherfile_name, p, plen);
	otherfile_name[plen] = '\0';
	if ((p = edf_findInHeader(header, "EDF_BinaryFilePosition")))
	    otherfile_skip = atoi(p);
    }

    if ((p = edf_findInHeader(header, "Dim_1")))
	dim1 = atoi(p);
    if ((p = edf_findInHeader(header, "Dim_2")))
	dim2 = atoi(p);
    if ((p = edf_findInHeader(header, "DataType"))) {
	k = lookup_table3_nth(edf_datatype_table, p);
	if (k < 0) { // unknown EDF DataType
	    error("Unknown EDF datatype \"%s\"", p);
	    free(header);
	    return retval;
	}
	datatype = edf_datatype_table[k].value;
	datalen = edf_datatype_table[k].sajzof;
    }
    int byteorder = THIS_ENDIAN_BYTE_ORDER;
    int swap_bytes = 0;
    if ((p = edf_findInHeader(header, "ByteOrder"))) {
	k = lookup_table_nth(edf_byteorder_table, p);
	if (k >= 0) {
	    byteorder = edf_byteorder_table[k].value;
	    swap_bytes = (byteorder != THIS_ENDIAN_BYTE_ORDER);
	}
    } else
	std::cerr <<"WARNING: ByteOrder not specified in the header! Not swapping bytes (figure may not be correct).\n";
    // Get and verify size of the data:
    int datasize = dim1 * dim2 * datalen;
    if ((p = edf_findInHeader(header, "Size"))) {
	int d = atoi(p);
	if (d != datasize) {
	    std::cerr <<"WARNING: Size " <<datasize <<" is not " <<dim1 <<'x' <<dim2
       		      <<"x" <<datalen <<" = " <<d <<". Supposing the latter.\n";
	}
    }
    int raster_axes = RASTER_AXES_XrightYdown; // default
    if ((p = edf_findInHeader(header, "RasterAxes"))) {
	raster_axes = lookup_table_nth(rasteraxes_table, p);
    }

#ifdef DEBUG
    std::cerr <<"Parsing header: dim1=" <<dim1 <<" dim2=" <<dim2
	      <<" headersize=" <<header_size <<" datasize=" <<datasize
	      <<"\n\t\tbyteorder=" <<byteorder <<" datatype=" <<datatype 
	      <<" datalen=" <<datalen <<" raster_axes =" <<raster_axes
	      <<"\n";
#endif

    // EHF files: binary data are in another file than the header file
    if (otherfile_name) {
	gzclose(inp);
	inp = gzopen(otherfile_name, "rb");
	if (!inp) {
	    error("Cannot open file \"%s\"", otherfile_name);
	    free(header);
	    return retval;
	}
	// gzrewind(inp);
	gzseek(inp, otherfile_skip, SEEK_SET);
    }

    retval(0) = octave_value(header);
    free(header);
    header = 0;

    if (nargout == 1)
	return retval;

    // read the data (image)

    // gzrewind(inp); gzread(inp, header, header_size); // ... not needed here

    char *image = new char [datasize];
    if (!image) {
	error("Not enough memory for the image\n");
	return retval;
    }

    if (datasize != gzread(inp, &image[0], datasize)) {
	error("The image cannot be read completely.");
	delete [] image;
	return retval;
    }

    if (swap_bytes) {
	char *pixel = image;
	for (int i=dim1*dim2; --i>=0; pixel += datalen)
	    byteswap( pixel, datalen );
    }
    gzclose(inp);

    // copy the image to octave's matrix object
    Matrix a(dim1, dim2);
    p = image;
    double d = 0;
    int yy = dim2-1;
    for (int y = 0; y < dim2; y++, yy--) {
    	for (int x = 0; x < dim1; x++) {
	    switch (datatype) {
#define X(T) d = *(T*)p; break;
		case U_CHAR_DATATYPE:	X(unsigned char);
		case CHAR_DATATYPE:	X(signed char);
		case U_SHORT_DATATYPE:	X(unsigned short);
		case SHORT_DATATYPE:	X(signed short);
		case U_INT_DATATYPE:	X(unsigned int);
		case INT_DATATYPE:	X(signed int);
		case U_L_INT_DATATYPE:	X(unsigned long int);
		case L_INT_DATATYPE:	X(signed long int);
		case FLOAT_DATATYPE:	X(float); 
		case DOUBLE_DATATYPE:	X(double);
#undef X
		case UNKNOWN_DATATYPE:	break; // make compiler happy
	    }
	    if (raster_axes == RASTER_AXES_XrightYup) 
		a(x,yy) = d; 
	    else // raster_axes == EDF_RASTER_AXES_XrightYdown
		a(x,y) = d; 
	    p += datalen;
	}
    }
    delete [] image;

    // return the image data
    if (nargout >= 2) {
	retval(1) = a;
    }
    return retval;
}

// eof pmedf_readC.cpp
