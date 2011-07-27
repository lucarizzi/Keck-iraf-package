#  Version Info: This file is distributed with version 2.401 of CFITSIO   */
# This file created from cfitsio/fitsio.h


# global variables */
 
define FLEN_FILENAME 1025 # max length of a filename  */
define FLEN_KEYWORD   72  # max length of a keyword (HIERARCH convention) */
define FLEN_CARD      81  # length of a FITS header card */
define FLEN_VALUE     71  # max length of a keyword value string */
define FLEN_COMMENT   73  # max length of a keyword comment string */
define FLEN_ERRMSG    81  # max length of a FITSIO error message */
define FLEN_STATUS    31  # max length of a FITSIO status text string */
 
define TBIT          1  # codes for FITS table data types */
define TBYTE        11
define TLOGICAL     14
define TSTRING      16
define TUSHORT      20
define TSHORT       21
define TUINT        30
define TINT         31
define TULONG       40
define TLONG        41
define TINT32BIT    41  # used when returning datatype of a column */
define TFLOAT       42
define TLONGLONG    81
define TDOUBLE      82
define TCOMPLEX     83
define TDBLCOMPLEX 163

define TYP_STRUC_KEY 10
define TYP_CMPRS_KEY 20
define TYP_SCAL_KEY  30
define TYP_NULL_KEY  40
define TYP_DIM_KEY   50
define TYP_RANG_KEY  60
define TYP_UNIT_KEY  70
define TYP_DISP_KEY  80
define TYP_HDUID_KEY 90
define TYP_CKSUM_KEY 100
define TYP_WCS_KEY   110
define TYP_REFSYS_KEY 120
define TYP_COMM_KEY  130
define TYP_CONT_KEY  140
define TYP_USER_KEY  150



define BYTE_IMG      8  # BITPIX code values for FITS image types */
define SHORT_IMG    16
define LONG_IMG     32
define LONGLONG_IMG 64
define FLOAT_IMG   -32
define DOUBLE_IMG  -64
                         # The following 2 codes are not true FITS         */
                         # datatypes; these codes are only used internally */
                         # within cfitsio to make it easier for users      */
                         # to deal with unsigned integers.                 */
define USHORT_IMG   20
define ULONG_IMG    40

define IMAGE_HDU  0  # Primary Array or IMAGE HDU */
define ASCII_TBL  1  # ASCII table HDU  */
define BINARY_TBL 2  # Binary table HDU */
define ANY_HDU   -1  # matches any HDU type */

define READONLY  0    # options when opening a file */
define READWRITE 1



#* Image compression algorithm types */
define MAX_COMPRESS_DIM     6
define RICE_1      11
define GZIP_1      21
define PLIO_1      31
define HCOMPRESS_1 41

#* error status codes */

define SKIP_TABLE       -104  ## move to 1st image when opening file */
define SKIP_IMAGE       -103  ## move to 1st table when opening file */
define SKIP_NULL_PRIMARY -102 ## skip null primary array when opening file */
define USE_MEM_BUFF     -101  ## use memory buffer when opening file */
define OVERFLOW_ERR      -11  ## overflow during datatype conversion */
define PREPEND_PRIMARY    -9  ## used in ffiimg to insert new primary array */
define SAME_FILE         101  ## input and output files are the same */
define TOO_MANY_FILES    103  ## tried to open too many FITS files */
define FILE_NOT_OPENED   104  ## could not open the named file */
define FILE_NOT_CREATED  105  ## could not create the named file */
define WRITE_ERROR       106  ## error writing to FITS file */
define END_OF_FILE       107  ## tried to move past end of file */
define READ_ERROR        108  ## error reading from FITS file */
define FILE_NOT_CLOSED   110  ## could not close the file */
define ARRAY_TOO_BIG     111  ## array dimensions exceed internal limit */
define READONLY_FILE     112  ## Cannot write to readonly file */
define MEMORY_ALLOCATION 113  ## Could not allocate memory */
define BAD_FILEPTR       114  ## invalid fitsfile pointer */
define NULL_INPUT_PTR    115  ## NULL input pointer to routine */
define SEEK_ERROR        116  ## error seeking position in file */

define BAD_URL_PREFIX    121  ## invalid URL prefix on file name */
define TOO_MANY_DRIVERS  122  ## tried to register too many IO drivers */
define DRIVER_INIT_FAILED 123  ## driver initialization failed */
define NO_MATCHING_DRIVER 124  ## matching driver is not registered */
define URL_PARSE_ERROR    125  ## failed to parse input file URL */
define RANGE_PARSE_ERROR  126  ## failed to parse input file URL */


define HEADER_NOT_EMPTY  201  ## header already contains keywords */
define KEY_NO_EXIST      202  ## keyword not found in header */
define KEY_OUT_BOUNDS    203  ## keyword record number is out of bounds */
define VALUE_UNDEFINED   204  ## keyword value field is blank */
define NO_QUOTE          205  ## string is missing the closing quote */
define BAD_KEYCHAR       207  ## illegal character in keyword name or card */
define BAD_ORDER         208  ## required keywords out of order */
define NOT_POS_INT       209  ## keyword value is not a positive integer */
define NO_END            210  ## couldn't find END keyword */
define BAD_BITPIX        211  ## illegal BITPIX keyword value*/
define BAD_NAXIS         212  ## illegal NAXIS keyword value */
define BAD_NAXES         213  ## illegal NAXISn keyword value */
define BAD_PCOUNT        214  ## illegal PCOUNT keyword value */
define BAD_GCOUNT        215  ## illegal GCOUNT keyword value */
define BAD_TFIELDS       216  ## illegal TFIELDS keyword value */
define NEG_WIDTH         217  ## negative table row size */
define NEG_ROWS          218  ## negative number of rows in table */
define COL_NOT_FOUND     219  ## column with this name not found in table */
define BAD_SIMPLE        220  ## illegal value of SIMPLE keyword  */
define NO_SIMPLE         221  ## Primary array doesn't start with SIMPLE */
define NO_BITPIX         222  ## Second keyword not BITPIX */
define NO_NAXIS          223  ## Third keyword not NAXIS */
define NO_NAXES          224  ## Couldn't find all the NAXISn keywords */
define NO_XTENSION       225  ## HDU doesn't start with XTENSION keyword */
define NOT_ATABLE        226  ## the CHDU is not an ASCII table extension */
define NOT_BTABLE        227  ## the CHDU is not a binary table extension */
define NO_PCOUNT         228  ## couldn't find PCOUNT keyword */
define NO_GCOUNT         229  ## couldn't find GCOUNT keyword */
define NO_TFIELDS        230  ## couldn't find TFIELDS keyword */
define NO_TBCOL          231  ## couldn't find TBCOLn keyword */
define NO_TFORM          232  ## couldn't find TFORMn keyword */
define NOT_IMAGE         233  ## the CHDU is not an IMAGE extension */
define BAD_TBCOL         234  ## TBCOLn keyword value < 0 or > rowlength */
define NOT_TABLE         235  ## the CHDU is not a table */
define COL_TOO_WIDE      236  ## column is too wide to fit in table */
define COL_NOT_UNIQUE    237  ## more than 1 column name matches template */
define BAD_ROW_WIDTH     241  ## sum of column widths not = NAXIS1 */
define UNKNOWN_EXT       251  ## unrecognizable FITS extension type */
define UNKNOWN_REC       252  ## unrecognizable FITS record */
define END_JUNK          253  ## END keyword is not blank */
define BAD_HEADER_FILL   254  ## Header fill area not blank */
define BAD_DATA_FILL     255  ## Data fill area not blank or zero */
define BAD_TFORM         261  ## illegal TFORM format code */
define BAD_TFORM_DTYPE   262  ## unrecognizable TFORM datatype code */
define BAD_TDIM          263  ## illegal TDIMn keyword value */
define BAD_HEAP_PTR      264  ## invalid BINTABLE heap address */

define BAD_HDU_NUM       301  ## HDU number < 1 or > MAXHDU */
define BAD_COL_NUM       302  ## column number < 1 or > tfields */
define NEG_FILE_POS      304  ## tried to move before beginning of file  */
define NEG_BYTES         306  ## tried to read or write negative bytes */
define BAD_ROW_NUM       307  ## illegal starting row number in table */
define BAD_ELEM_NUM      308  ## illegal starting element number in vector */
define NOT_ASCII_COL     309  ## this is not an ASCII string column */
define NOT_LOGICAL_COL   310  ## this is not a logical datatype column */
define BAD_ATABLE_FORMAT 311  ## ASCII table column has wrong format */
define BAD_BTABLE_FORMAT 312  ## Binary table column has wrong format */
define NO_NULL           314  ## null value has not been defined */
define NOT_VARI_LEN      317  ## this is not a variable length column */
define BAD_DIMEN         320  ## illegal number of dimensions in array */
define BAD_PIX_NUM       321  ## first pixel number greater than last pixel */
define ZERO_SCALE        322  ## illegal BSCALE or TSCALn keyword = 0 */
define NEG_AXIS          323  ## illegal axis length < 1 */

#define NOT_GROUP_TABLE         340
#define HDU_ALREADY_MEMBER      341
#define MEMBER_NOT_FOUND        342
#define GROUP_NOT_FOUND         343
#define BAD_GROUP_ID            344
#define TOO_MANY_HDUS_TRACKED   345
#define HDU_ALREADY_TRACKED     346
#define BAD_OPTION              347
#define IDENTICAL_POINTERS      348
#define BAD_GROUP_ATTACH        349
#define BAD_GROUP_DETACH        350

#define BAD_I2C           401  ## bad int to formatted string conversion */
#define BAD_F2C           402  ## bad float to formatted string conversion */
#define BAD_INTKEY        403  ## can't interprete keyword value as integer */
#define BAD_LOGICALKEY    404  ## can't interprete keyword value as logical */
#define BAD_FLOATKEY      405  ## can't interprete keyword value as float */
#define BAD_DOUBLEKEY     406  ## can't interprete keyword value as double */
#define BAD_C2I           407  ## bad formatted string to int conversion */
#define BAD_C2F           408  ## bad formatted string to float conversion */
#define BAD_C2D           409  ## bad formatted string to double conversion */
#define BAD_DATATYPE      410  ## bad keyword datatype code */
#define BAD_DECIM         411  ## bad number of decimal places specified */
#define NUM_OVERFLOW      412  ## overflow during datatype conversion */

# define DATA_COMPRESSION_ERR 413  ## error in imcompress routines */
# define DATA_DECOMPRESSION_ERR 414 ## error in imcompress routines */
# define NO_COMPRESSED_TILE  415 ## compressed tile doesn't exist */

#define BAD_DATE          420  ## error in date or time conversion */

#define PARSE_SYNTAX_ERR  431  ## syntax error in parser expression */
#define PARSE_BAD_TYPE    432  ## expression did not evaluate to desired type */
#define PARSE_LRG_VECTOR  433  ## vector result too large to return in array */
#define PARSE_NO_OUTPUT   434  ## data parser failed not sent an out column */
#define PARSE_BAD_COL     435  ## bad data encounter while parsing column */
#define PARSE_BAD_OUTPUT  436  ## Output file not of proper type          */

#define ANGLE_TOO_BIG     501  ## celestial angle too large for projection */
#define BAD_WCS_VAL       502  ## bad celestial coordinate or pixel value */
#define WCS_ERROR         503  ## error in celestial coordinate calculation */
#define BAD_WCS_PROJ      504  ## unsupported type of celestial projection */
#define NO_WCS_KEY        505  ## celestial coordinate keywords not found */
#define APPROX_WCS_KEY    506  ## approximate WCS keywords were calculated */

#define NO_CLOSE_ERROR    999  ## special value used internally to switch off */
                               ## the error message from ffclos and ffchdu */

