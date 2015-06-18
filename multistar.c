

#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PI 3.14159265358979323846

#include "uthash.h"

#include "postgres.h"
#include "funcapi.h"
#include "fmgr.h"

#include "access/htup_details.h" // NO TUP WARNINGS
#include "libpq/pqformat.h"
#include "catalog/pg_type.h"
#include "executor/spi.h"
#include "executor/executor.h"  /* for GetAttributeByName() */
#include "utils/array.h"
#include "utils/lsyscache.h"


#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif


/* Triangle contains 3 ids to its 3 vertices */
typedef struct Triangle {
  int v0;
  int v1;
  int v2;

} Triangle;

/* Vector used for normals, slopes etc. */
typedef struct Vec {
  float i;
  float j;
  float k;
} Vec;

/* Point struct */
typedef struct
{
    float4  x, y, z;
} Point;

/* Taken from PostGIS BOX3D to parse bounding box 
https://github.com/postgis/postgis/blob/svn-trunk/liblwgeom/liblwgeom.h.in
line 280 */ 
typedef struct BOX3D {
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  int32_t srid;
} BOX3D;
  
/* Vstar is a struct used to store the values from the fetched rows */
typedef struct Vstar {
  int id;
  int patchid;
  int sid;            // stitched id from pairing id & patchid
  float8 coords[3];
  int *star;
  int *dim;
  bool flag;
  UT_hash_handle hh;  // makes this structure hashable
} Vstar;

/* Hash table of Vstars */
struct Vstar *htable = NULL;

/* Multistar structure, stores complete row from db */
typedef struct Multistar {
  int patchid;
  int nstars;
  int32 *offsets;
  float4 *points;
  int32 *stars;
  BOX3D *bbox;
} Multistar;

/* Stack structure */
typedef struct stack {
  int data;
  struct stack *next;
} STACK;

/* Structure for packing and unpacking with pairing function */
typedef struct Link {
  uint16_t patch;
  uint16_t point;
} Link;

/*#define SRF_RETURN_NEXT_NULL(_funcctx)
           do {
                   ReturnSetInfo      *rsi;
                   (_funcctx)->call_cntr++;
                   rsi = (ReturnSetInfo *) fcinfo->resultinfo;
                   rsi->isDone = ExprMultipleResult;
                   PG_RETURN_NULL();
           } while (0) */





typedef struct
{
    double x, y, z;
} Pointd;

typedef struct
{
  int32 a, b, c, n1, n2, n3;
} TTriangle;


Datum npoints(PG_FUNCTION_ARGS);
Datum getx(PG_FUNCTION_ARGS);
Datum getstar(PG_FUNCTION_ARGS);
Datum getpoint(PG_FUNCTION_ARGS);
Datum getpoints(PG_FUNCTION_ARGS);
Datum gettin(PG_FUNCTION_ARGS);
Datum gettin_convex(PG_FUNCTION_ARGS);
Datum getrow(PG_FUNCTION_ARGS);

Datum simplify(PG_FUNCTION_ARGS);

char machine_endian(void);


 


///-- added functions for the Triangle type
Datum triangle_in(PG_FUNCTION_ARGS);
Datum triangle_out(PG_FUNCTION_ARGS);
Datum triangle_recv(PG_FUNCTION_ARGS);
Datum triangle_send(PG_FUNCTION_ARGS);



Datum point_location(PG_FUNCTION_ARGS);
Datum pl(PG_FUNCTION_ARGS);
Datum range_query(PG_FUNCTION_ARGS);
Datum estimate_z(PG_FUNCTION_ARGS);
Datum profile(PG_FUNCTION_ARGS);
Datum profile_count_intersections(PG_FUNCTION_ARGS);
Datum testmulti(PG_FUNCTION_ARGS);  

int get_starting_id(float8 dest[], int patchid);
int get_starting_id_convex(float8 dest[], int patchid);

Multistar* fetch_row(void* pplan2, int patchid);
Vstar* fetch_row_vstar(int id, void* pplan, int patchid);
Vstar* fetch_row_vstar2(int id, Multistar **row, void* pplan_row, int patchid);

void* trianglez_to_wkb(Vstar **tr, size_t *wkbsize);
Vstar** get_starting_triangle_convex(float8 a[], float8 b[], int patchid);
Vstar** jump_and_walk(float8 dest[], void *pplan2, int startid, float8 power, int patchid);
void walk_straight(Vstar **tr, Multistar **multistar, float8 a[], float8 b[], float8 bboxmin[], float8 bboxmax[], void *pplan, int patchid);
void* walk_straight_and_report_intersections(Vstar **tr, float8 a[], float8 b[], void *pplan, bool justcount, int patchid, size_t *size);
int orient2d(float8 *a, float8 *b, float8 *c);
bool intersect(float8 *a, float8 *b, float8 *c, float8 *d);
void get_intersection_z(float8 *a, float8 *b, float8 *c, float8 *d, double *x);
double distance_sqr(float8 *a, float8 *b);
int find_previous_in_star(int star[], int dim, int v);
int find_next_in_star(int star[], int dim, int v);
void* prepare_plan(void);
void* prepare_plan2(void);
int no_points(int patchid); 
double area_triangle(float8 *a, float8 *b, float8 *c);
bool inside_rectangle(float8 *p, float8 *bboxmin, float8 *bboxmax);
void push(STACK **head, int value);
int pop(STACK **head);
uint32_t stitch(uint16_t patch_id, uint16_t point_id);
Link* unstitch(uint32_t packed);

Vec* normal(Point *p[]);
float aspect(Vec *normal);
float slope(Vec *normal);


uint32_t stitch(uint16_t patch_id, uint16_t point_id)
{
  uint32_t *packed = palloc0(sizeof(uint32_t));
  *packed = (uint32_t) ((point_id << 16) | patch_id);
  ////elog(INFO,"Stuff in %i %i, out %i, ",patch_id,point_id,packed);
  return *packed;
}

Link* unstitch(uint32_t packed)
{
  Link *link;
  link = palloc(sizeof(Link));

  ////elog(INFO,"Packed %i",packed);
  link->point = (packed >> 16);
  link->patch = (packed & 0xFFFF);
  ////elog(INFO,"Unpacked %i %i",link->patch,link->point);
  return link;
}

char machine_endian(void)
{
  static int check_int = 1; /* dont modify this!!! */
  return *((char *) &check_int); /* 0 = big endian | xdr,
                                 * 1 = little endian | ndr
                                   */
}


// FUNCTION TO GET NUMBER OF POINTS 
// INPUT TAKES (TREE)
PG_FUNCTION_INFO_V1(npoints);
Datum
npoints(PG_FUNCTION_ARGS)
{
  ArrayType  *input;
    int32 *n;
    input = PG_GETARG_ARRAYTYPE_P(0);
  n = (int32*) ARR_DATA_PTR(input);
    PG_RETURN_INT32(abs(n[0]));
}

// GET POINT FUNCTION
// INPUT TAKES (i, TREE, BYTEA)
PG_FUNCTION_INFO_V1(getx);
Datum 
getx(PG_FUNCTION_ARGS)
{
  ArrayType *input;
  int32 *n;
  int32 len;
  Point *points;
  int32 i;
  float8 x;

  // get number of points
  // to slice bytea
  input = PG_GETARG_ARRAYTYPE_P(1);
  n = (int32*) ARR_DATA_PTR(input);
  len = abs(n[0]);

  // PROBABLY NEEDS VARLEN HEADER TO ADD
  i = PG_GETARG_INT32(0);
  if (!(i > 0 && i < len)) {
    PG_RETURN_NULL();
  }
  else {
  points = (Point*) PG_GETARG_BYTEA_P_SLICE(2,VARHDRSZ,len*sizeof(Point)+VARHDRSZ);
    x = points[i].x;

    PG_RETURN_FLOAT8(x);
    }
}

PG_FUNCTION_INFO_V1(getpoint);
Datum 
getpoint(PG_FUNCTION_ARGS)
{
  //ArrayType *input;
  //int32 *n;
  //int32 len;
  Point *points;
  int32 i;
  float8 x;
  float8 y;
  float8 z;

    TupleDesc resultTupleDesc; 
    Oid resultTypeId; 
    Datum retvals[3]; 
    bool  retnulls[3]; 
    HeapTuple rettuple; 

  // get number of points
  // to slice bytea
  //input = PG_GETARG_ARRAYTYPE_P(1);
  //n = (int32*) ARR_DATA_PTR(input);
  //len = abs(n[0]);

  // PROBABLY NEEDS VARLEN HEADER TO ADD
  i = PG_GETARG_INT32(0)-1;
  //if (!(i > 0 && i < len)) {
  //  PG_RETURN_NULL();
  //}
  //else {
  //points = (Point*) PG_GETARG_BYTEA_P(2);
  points = (Point*) ( ( (uint8_t*) PG_GETARG_BYTEA_P(2) ) + 4);
    x = points[i].x;
    y = points[i].y;
    z = points[i].z;


    ereport(INFO, 
               (errmsg("arg: (x: %f,y: %f, z: %f var: %i len: %li)", x, y, z, VARHDRSZ, sizeof(Point))) ); 


    get_call_result_type(fcinfo, &resultTypeId, &resultTupleDesc); 
    Assert(resultTypeId == TYPEFUNC_COMPOSITE); 
    BlessTupleDesc(resultTupleDesc);  

    retvals[0] = Float8GetDatum(x);
    retvals[1] = Float8GetDatum(y);
    retvals[2] = Float8GetDatum(z);

    retnulls[0] = false;
    retnulls[1] = false;
    retnulls[2] = false;

    rettuple = heap_form_tuple( resultTupleDesc, retvals, retnulls ); 

    PG_RETURN_DATUM( HeapTupleGetDatum( rettuple ) );
}

PG_FUNCTION_INFO_V1(getstar);
Datum 
getstar(PG_FUNCTION_ARGS)
{
  ArrayType *input, *result;
  Datum *results;
  int32 *n;
  int32 len, starlen;
  int32 *stars;
  int32 i, a, b;
  Oid         eltype;
    int16       elmlen;
    bool        elmbyval;
    char        elmalign;

  // get number of points
  // to slice bytea
  input = PG_GETARG_ARRAYTYPE_P(1);
  n = (int32*) ARR_DATA_PTR(input);
  len = abs(n[0]);

  // PROBABLY NEEDS VARLEN HEADER TO ADD
  i = PG_GETARG_INT32(0);

  if (!(i > 0 && i <= len)) {
    PG_RETURN_NULL();
  }
  // END OF LIST
  else if (i == len) {
    b = VARSIZE(PG_GETARG_BYTEA_P(2))-VARHDRSZ-len*sizeof(Point);
    a = abs(n[i]);
    starlen = (abs(b)-abs(a))/sizeof(int32);
  }
  // USING OFFSET LIST
  else {
    a = abs(n[i]);
    b = abs(n[i+1]);
    starlen = (abs(b)-abs(a))/sizeof(int32);
  }
  
  ereport(INFO, 
            (errmsg("returning %i items in star\n", starlen)));

  stars = (int32*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,len*sizeof(Point)+a,starlen*sizeof(int32))) + 4);

    
    eltype = ARR_ELEMTYPE(input);
    get_typlenbyvalalign(eltype, &elmlen, &elmbyval, &elmalign);
  results = (Datum *)palloc(starlen * sizeof(Datum));

  for (i=0;i<starlen;i++) {
    results[i] = DatumGetInt32(stars[i]);
  }

  result = construct_array(results, starlen,
                                eltype, elmlen, elmbyval, elmalign);
  PG_RETURN_ARRAYTYPE_P(result);

}

PG_FUNCTION_INFO_V1(getrow);
Datum 
getrow(PG_FUNCTION_ARGS)
{
  ArrayType *input, *result;
  Datum *results;
  int32 *n;
  int32 len, starlen;
  int32 *stars;
  int32 i, l, a, b;
  Oid         eltype;
    int16       elmlen;
    bool        elmbyval;
    char        elmalign;
    Point *points;
  float8 x;
  float8 y;
  float8 z;
  TupleDesc resultTupleDesc; 
    Oid resultTypeId; 
    Datum retvals[5]; 
    bool  retnulls[5]; 
    HeapTuple rettuple; 

  // get number of points
  // to slice bytea
  input = PG_GETARG_ARRAYTYPE_P(1);
  n = (int32*) ARR_DATA_PTR(input);
  len = abs(n[0]);
  //elog(INFO,"GR Wanted is %i",i);

  // PROBABLY NEEDS VARLEN HEADER TO ADD
  i = PG_GETARG_INT32(0);

  if (!(i > 0 && i <= len)) {
    PG_RETURN_NULL();
  }
  // END OF LIST
  else if (i == len) {
    b = VARSIZE(PG_GETARG_BYTEA_P(3))-VARHDRSZ;
    a = abs(n[i]);
    starlen = (abs(b)-abs(a))/sizeof(int32);
  }
  // USING OFFSET LIST
  else {
    a = abs(n[i]);
    b = abs(n[i+1]);
    starlen = (abs(b)-abs(a))/sizeof(int32);
  }
  
  stars = (int32*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(3,a,starlen*sizeof(int32)) + 4));
  //points = (Point*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,0,VARSIZE(PG_GETARG_BYTEA_P(2)))) +4);
  points = (Point*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,(i-1)*sizeof(Point),i*sizeof(Point))) +4);

    x = points[0].x;
    y = points[0].y;
    z = points[0].z;

    eltype = ARR_ELEMTYPE(input);
    get_typlenbyvalalign(eltype, &elmlen, &elmbyval, &elmalign);
  results = (Datum *)palloc(starlen * sizeof(Datum));

  for (l=0;l<starlen;l++) {
    results[l] = DatumGetInt32(stars[l]);
  }

  result = construct_array(results, starlen,
                                eltype, elmlen, elmbyval, elmalign);

  get_call_result_type(fcinfo, &resultTypeId, &resultTupleDesc); 
    Assert(resultTypeId == TYPEFUNC_COMPOSITE); 
    BlessTupleDesc(resultTupleDesc);  

    retvals[0] = Int32GetDatum(i);
    retvals[1] = Float8GetDatum(x);
    retvals[2] = Float8GetDatum(y);
    retvals[3] = Float8GetDatum(z);
    retvals[4] = PointerGetDatum(result);

    retnulls[0] = false;
    retnulls[1] = false;
    retnulls[2] = false;
    retnulls[3] = false;
    retnulls[4] = false;

    rettuple = heap_form_tuple( resultTupleDesc, retvals, retnulls ); 
    PG_RETURN_DATUM( HeapTupleGetDatum( rettuple ) );
}


PG_FUNCTION_INFO_V1(getpoints);
Datum 
getpoints(PG_FUNCTION_ARGS)
{
  ArrayType *input;
  FuncCallContext     *funcctx;
  int32 *n;
  int32 len;
  int32 call_cntr=0;
  Point *points;

  if (SRF_IS_FIRSTCALL())
    {

        MemoryContext   oldcontext;
    funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

    input = PG_GETARG_ARRAYTYPE_P(0);
    n = (int32*) ARR_DATA_PTR(input);
    len = abs(n[0]);
    funcctx->max_calls = len;

    //points = palloc0(sizeof(points));
    //points = (Point*) PG_GETARG_BYTEA_P_SLICE(1,0,len*sizeof(Point));
    //points = (Point*) ( ( (uint8_t*) PG_GETARG_BYTEA_PP(2) ) + 4);
    points = (Point*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(1,0,len*sizeof(Point)) ) + 4);

    funcctx->user_fctx = points;

        if (get_call_result_type(fcinfo, NULL, &funcctx->tuple_desc) != TYPEFUNC_COMPOSITE)
            ereport(ERROR,
                    (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                     errmsg("function returning record called in context "
                            "that cannot accept type record")));
      BlessTupleDesc(funcctx->tuple_desc);  
    MemoryContextSwitchTo(oldcontext);

  }

    funcctx = SRF_PERCALL_SETUP();
    points = funcctx->user_fctx;
    call_cntr = funcctx->call_cntr;
 
  if (call_cntr < funcctx->max_calls) {

    TupleDesc resultTupleDesc; 
      Oid resultTypeId; 
      Datum retvals[4]; 
      bool  retnulls[4]; 
      HeapTuple rettuple;
      get_call_result_type(fcinfo, &resultTypeId, &resultTupleDesc); 
      Assert(resultTypeId == TYPEFUNC_COMPOSITE); 
      retvals[0] = Int32GetDatum(call_cntr+1);
      retvals[1] = Float8GetDatum(points[call_cntr].x);
      retvals[2] = Float8GetDatum(points[call_cntr].y);
      retvals[3] = Float8GetDatum(points[call_cntr].z);
      retnulls[0] = false;
      retnulls[1] = false;
      retnulls[2] = false;
      retnulls[3] = false;
      rettuple = heap_form_tuple( funcctx->tuple_desc, retvals, retnulls ); 
      SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum( rettuple ));
    }
    else {
      SRF_RETURN_DONE(funcctx);
    }
}


typedef struct { 
  ArrayType *offset; 
  Point *points; 
  int32 *stars; 
  } my_user_fctx;

PG_FUNCTION_INFO_V1(gettin);
Datum 
gettin(PG_FUNCTION_ARGS)
{
  FuncCallContext     *funcctx;
  int32 *n,*stars;
  int32 call_cntr=0;
  Point *points;
  ArrayType *input, *result;
  Datum *results;
  int32 len, starlen;
  my_user_fctx *ref = palloc(sizeof(my_user_fctx));
  Oid         eltype;
    int16       elmlen;
    bool        elmbyval;
    char        elmalign;

  if (SRF_IS_FIRSTCALL())
    {

        MemoryContext   oldcontext;
    funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

    input = PG_GETARG_ARRAYTYPE_P(0);
    //n = palloc(VARSIZE(input)-VARHDRSZ);
    n = (int32*) ARR_DATA_PTR(input);
    len = abs(n[0]);

    //points = palloc(len*sizeof(Point));
    points = (Point*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(1,0,len*sizeof(Point)+4 )) + 4);
    stars = palloc(VARSIZE(PG_GETARG_BYTEA_P(2))-VARHDRSZ);
    stars = (int32*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,0,VARSIZE(PG_GETARG_BYTEA_P(2)))) +4);
    funcctx->max_calls = len;

    ref->offset = input;
    ref->points = points;
    ref->stars = stars;
    funcctx->user_fctx = ref;

        if (get_call_result_type(fcinfo, NULL, &funcctx->tuple_desc) != TYPEFUNC_COMPOSITE)
            ereport(ERROR,
                    (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                     errmsg("function returning record called in context "
                            "that cannot accept type record")));
      BlessTupleDesc(funcctx->tuple_desc);  
    MemoryContextSwitchTo(oldcontext);
  }

    funcctx = SRF_PERCALL_SETUP();
    call_cntr = funcctx->call_cntr;

  if (call_cntr < funcctx->max_calls) {
    TupleDesc resultTupleDesc; 
      Oid resultTypeId; 
      Datum retvals[5]; 
      bool  retnulls[5]; 
      HeapTuple rettuple;
      int32 i, a, b;
    my_user_fctx *ref = palloc(sizeof(my_user_fctx));

    ref = funcctx->user_fctx;

      get_call_result_type(fcinfo, &resultTypeId, &resultTupleDesc); 
      Assert(resultTypeId == TYPEFUNC_COMPOSITE); 

    n = (int32*) ARR_DATA_PTR(ref->offset);
    points = (Point*) ref->points;
    stars = (int32*) ref->stars;

    if (call_cntr == funcctx->max_calls-1) {
      b = (VARSIZE(PG_GETARG_BYTEA_P(2))-VARHDRSZ)/sizeof(int32);
      a = abs(n[call_cntr+1])/sizeof(int32);
      starlen = (abs(b)-abs(a));
    }
    // USING OFFSET LIST
    else {
      a = abs(n[call_cntr+1])/sizeof(int32);
      b = abs(n[call_cntr+2])/sizeof(int32);
      starlen = (abs(b)-abs(a));
    }
    // I KNOW, A and B are wrong!
    //ereport(INFO, 
      //      (errmsg("a,b items %i %i %i\n", a,b,starlen)));
    //stars = (int32*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,a,starlen*sizeof(int32))) + 4);
   
      eltype = ARR_ELEMTYPE(ref->offset);
      get_typlenbyvalalign(eltype, &elmlen, &elmbyval, &elmalign);

    results = (Datum *)palloc(starlen * sizeof(Datum));
    for (i=0;i<starlen;i++) {
      results[i] = DatumGetInt32(stars[a+i]);
    }
    result = construct_array(results, starlen,
                                  eltype, elmlen, elmbyval, elmalign);

      retvals[0] = Int32GetDatum(call_cntr+1);
      retvals[1] = Float8GetDatum((float8)points[call_cntr].x);
      retvals[2] = Float8GetDatum((float8)points[call_cntr].y);
      retvals[3] = Float8GetDatum((float8)points[call_cntr].z);
      retvals[4] = PointerGetDatum(result);

      retnulls[0] = false;
      retnulls[1] = false;
      retnulls[2] = false;
      retnulls[3] = false;
      retnulls[4] = false;
      rettuple = heap_form_tuple( funcctx->tuple_desc, retvals, retnulls ); 
      SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum( rettuple ));
    }
    else {
      SRF_RETURN_DONE(funcctx);
    }
}

// GET CURRENT QUERY + LOCAL ID
// REPLACE SQL QUERY WITH gettin() WHERE ID IS
// NEXT TIME IN LOCAL PATCHES, RUN LOCAL GETTIN AND WALK AROUND IN THERE
// WITHOUT NEEDING SPI

// UPGRADE THIS FUNCTION TO ACTUALLY USE THE 9.4 RETURN EMPTY CALL
PG_FUNCTION_INFO_V1(gettin_convex);
Datum 
gettin_convex(PG_FUNCTION_ARGS)
{
  FuncCallContext     *funcctx;
  int32 *n,*stars;
  int32 call_cntr=0;
  Point *points;
  ArrayType *input, *result;
  Datum *results;
  int32 len, starlen;
  my_user_fctx *ref = palloc(sizeof(my_user_fctx));
  Oid         eltype;
    int16       elmlen;
    bool        elmbyval;
    char        elmalign;
    int convex;

  if (SRF_IS_FIRSTCALL())
    {

        MemoryContext   oldcontext;
    funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

    input = PG_GETARG_ARRAYTYPE_P(0);
    //n = palloc(VARSIZE(input)-VARHDRSZ);
    n = (int32*) ARR_DATA_PTR(input);
    len = abs(n[0]);

    //points = palloc(len*sizeof(Point));
    points = (Point*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(1,0,len*sizeof(Point)+4 )) + 4);
    stars = palloc(VARSIZE(PG_GETARG_BYTEA_P(2))-VARHDRSZ);
    stars = (int32*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,0,VARSIZE(PG_GETARG_BYTEA_P(2)))) +4);
    funcctx->max_calls = len;

    ref->offset = input;
    ref->points = points;
    ref->stars = stars;
    funcctx->user_fctx = ref;

        if (get_call_result_type(fcinfo, NULL, &funcctx->tuple_desc) != TYPEFUNC_COMPOSITE)
            ereport(ERROR,
                    (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                     errmsg("function returning record called in context "
                            "that cannot accept type record")));
      BlessTupleDesc(funcctx->tuple_desc);  
    MemoryContextSwitchTo(oldcontext);
  }

    funcctx = SRF_PERCALL_SETUP();
    call_cntr = funcctx->call_cntr;

    //while (convex > 0)
    //{
    //  funcctx->call_cntr++;
    //  //elog(INFO,"Skipping this non convex id %i",convex);
    //  convex = (n[funcctx->call_cntr+1])/sizeof(int32);
    //}

  if (funcctx->call_cntr < funcctx->max_calls) {
    TupleDesc resultTupleDesc; 
      Oid resultTypeId; 
      Datum retvals[5]; 
      bool  retnulls[5]; 
      HeapTuple rettuple;
      int32 i, a, b;
      //int32 convex;

    my_user_fctx *ref = palloc(sizeof(my_user_fctx));
      ref = funcctx->user_fctx;
      n = (int32*) ARR_DATA_PTR(ref->offset);
      //convex = (n[call_cntr+1])/sizeof(int32);
      //elog(INFO,"Convex now %i",convex);

      get_call_result_type(fcinfo, &resultTypeId, &resultTupleDesc); 
      Assert(resultTypeId == TYPEFUNC_COMPOSITE); 

    
    points = (Point*) ref->points;
    stars = (int32*) ref->stars;

    if (call_cntr == funcctx->max_calls-1) {
      convex = (n[call_cntr+1])/sizeof(int32);
      b = (VARSIZE(PG_GETARG_BYTEA_P(2))-VARHDRSZ)/sizeof(int32);
      a = abs(n[call_cntr+1])/sizeof(int32);
      starlen = (abs(b)-abs(a));
    }
    // USING OFFSET LIST
    else {
      a = abs(n[call_cntr+1])/sizeof(int32);
      b = abs(n[call_cntr+2])/sizeof(int32);
      starlen = (abs(b)-abs(a));
      convex = (n[call_cntr+1])/sizeof(int32);
    }
    // I KNOW, A and B are wrong!
    //ereport(INFO, 
      //      (errmsg("a,b items %i %i %i\n", a,b,starlen)));
    //stars = (int32*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,a,starlen*sizeof(int32))) + 4);
   
      if (convex < 0)
      {
        eltype = ARR_ELEMTYPE(ref->offset);
        get_typlenbyvalalign(eltype, &elmlen, &elmbyval, &elmalign);

      results = (Datum *)palloc(starlen * sizeof(Datum));
      for (i=0;i<starlen;i++) {
        results[i] = DatumGetInt32(stars[a+i]);
      }
      result = construct_array(results, starlen,
                                    eltype, elmlen, elmbyval, elmalign);

        retvals[0] = Int32GetDatum(call_cntr+1);
        retvals[1] = Float8GetDatum((float8)points[call_cntr].x);
        retvals[2] = Float8GetDatum((float8)points[call_cntr].y);
        retvals[3] = Float8GetDatum((float8)points[call_cntr].z);
        retvals[4] = PointerGetDatum(result);

        retnulls[0] = false;
        retnulls[1] = false;
        retnulls[2] = false;
        retnulls[3] = false;
        retnulls[4] = false;
        rettuple = heap_form_tuple( funcctx->tuple_desc, retvals, retnulls ); 
        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum( rettuple ));
      }
      else
      {
        retnulls[0] = true;
        retnulls[1] = true;
        retnulls[2] = true;
        retnulls[3] = true;
      retnulls[4] = true;
        rettuple = heap_form_tuple( funcctx->tuple_desc, retvals, retnulls ); 
        SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum( rettuple ));
        //SRF_RETURN_NEXT_NULL(funcctx);
      }

    }
    else {
      SRF_RETURN_DONE(funcctx);
    }
}

PG_FUNCTION_INFO_V1(simplify);
Datum 
simplify(PG_FUNCTION_ARGS)
{
  FuncCallContext     *funcctx;
  int32 *n,*stars;
  int32 call_cntr=0;
  Point *points;
  ArrayType *input, *result;
  Datum *results;
  int32 len, starlen;
  my_user_fctx *ref = palloc(sizeof(my_user_fctx));
  Oid         eltype;
    int16       elmlen;
    bool        elmbyval;
    char        elmalign;
    float dif;

  if (SRF_IS_FIRSTCALL())
    {

        MemoryContext   oldcontext;
    funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

    input = PG_GETARG_ARRAYTYPE_P(0);
    //n = palloc(VARSIZE(input)-VARHDRSZ);
    n = (int32*) ARR_DATA_PTR(input);
    len = abs(n[0]);

    //points = palloc(len*sizeof(Point));
    points = (Point*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(1,0,len*sizeof(Point)+4 )) + 4);
    stars = palloc(VARSIZE(PG_GETARG_BYTEA_P(2))-VARHDRSZ);
    stars = (int32*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,0,VARSIZE(PG_GETARG_BYTEA_P(2)))) +4);
    funcctx->max_calls = len;

    ref->offset = input;
    ref->points = points;
    ref->stars = stars;
    funcctx->user_fctx = ref;

        if (get_call_result_type(fcinfo, NULL, &funcctx->tuple_desc) != TYPEFUNC_COMPOSITE)
            ereport(ERROR,
                    (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                     errmsg("function returning record called in context "
                            "that cannot accept type record")));
      BlessTupleDesc(funcctx->tuple_desc);  
    MemoryContextSwitchTo(oldcontext);
  }

    funcctx = SRF_PERCALL_SETUP();
    call_cntr = funcctx->call_cntr;

  if (call_cntr < funcctx->max_calls) {
    TupleDesc resultTupleDesc; 
      Oid resultTypeId; 
      Datum retvals[4]; 
      bool  retnulls[4]; 
      HeapTuple rettuple;
      int32 i, a, b;
    my_user_fctx *ref = palloc(sizeof(my_user_fctx));

    ref = funcctx->user_fctx;

      get_call_result_type(fcinfo, &resultTypeId, &resultTupleDesc); 
      Assert(resultTypeId == TYPEFUNC_COMPOSITE); 

    n = (int32*) ARR_DATA_PTR(ref->offset);
    points = (Point*) ref->points;
    stars = (int32*) ref->stars;

    if (call_cntr == funcctx->max_calls-1) {
      b = (VARSIZE(PG_GETARG_BYTEA_P(2))-VARHDRSZ)/sizeof(int32);
      a = abs(n[call_cntr+1])/sizeof(int32);
      starlen = (abs(b)-abs(a));
    }
    // USING OFFSET LIST
    else {
      a = abs(n[call_cntr+1])/sizeof(int32);
      b = abs(n[call_cntr+2])/sizeof(int32);
      starlen = (abs(b)-abs(a));
    }
    // I KNOW, A and B are wrong!
    //ereport(INFO, 
      //      (errmsg("a,b items %i %i %i\n", a,b,starlen)));
    //stars = (int32*) ( ( (uint8_t*) PG_GETARG_BYTEA_P_SLICE(2,a,starlen*sizeof(int32))) + 4);
   
      /*eltype = ARR_ELEMTYPE(ref->offset);
      get_typlenbyvalalign(eltype, &elmlen, &elmbyval, &elmalign);

    results = (Datum *)palloc(starlen * sizeof(Datum));
    for (i=0;i<starlen;i++) {
      results[i] = DatumGetInt32(stars[a+i]);
    }
    result = construct_array(results, starlen,
                                  eltype, elmlen, elmbyval, elmalign);
    */

    dif = 0; 
    for (i=0;i<starlen;i++) {
      /* Always return convex hull! */
      if (stars[a+i] == 0) {
        dif = 999;
        break;
      }
      /* Ignore bucket edges for now */
      else if (stars[a+i] > 65000) {
        dif = 0;
        break;
      }
      else {
        dif += ((float8)points[call_cntr].z - (float8) points[stars[a+i]].z);
      }
    }   

      //retvals[0] = Int32GetDatum(call_cntr+1);
      retvals[0] = Float8GetDatum((float8)points[call_cntr].x);
      retvals[1] = Float8GetDatum((float8)points[call_cntr].y);
      retvals[2] = Float8GetDatum((float8)points[call_cntr].z);
      retvals[3] = Float8GetDatum(dif);

      retnulls[0] = false;
      retnulls[1] = false;
      retnulls[2] = false;
      retnulls[3] = false;
      //retnulls[4] = false;
      rettuple = heap_form_tuple( funcctx->tuple_desc, retvals, retnulls ); 
      SRF_RETURN_NEXT(funcctx, HeapTupleGetDatum( rettuple ));
    }
    else {
      SRF_RETURN_DONE(funcctx);
    }
}

/* ******************************************************* */
/// CREATE A STACK OR BFS
void push(STACK **head, int value)
{
  STACK *node = (STACK*) malloc(sizeof(STACK));  /* create a new node */
  node->data = value;
  if (*head == NULL) { 
    node->next = NULL;
  }
  else 
    node->next = *head;
  *head = node;
}

int pop(STACK **head)
{
  if (*head == NULL) {                          /* stack is empty */
    return -1;
  } 
  else {                                     /* pop a node */
    STACK *top = *head;
    int value = top->data;
    *head = top->next;
    free(top);
    return value;
  }
}

/// RETURNS DISTANCE SQUARED
/// SAVE ROOTING SOMETHING
double distance_sqr(float8 *a, float8 *b)
{
  ////elog(INFO, "Distance of points a: %f %f and b: %f %f ",a[0],a[1],b[0],b[1]);
  return ( (a[0] - b[0]) * (a[0] - b[0]) ) + ( (a[1] - b[1]) * (a[1] - b[1]) );
}

int no_points(int patchid) {
  char sql[100];
  HeapTuple row;
  TupleDesc rdesc;
  int32 re;
  bool isnull = false;

  sprintf(sql, "SELECT @tree[1] from multistar_g37_l8 WHERE id = %i",patchid);
  SPI_execute(sql, true, 1);
  rdesc = SPI_tuptable->tupdesc;
  row = SPI_tuptable->vals[0];
  re = DatumGetInt32(SPI_getbinval(row, rdesc, 1, &isnull));
  return (int) re;
}

/// PREPARE SPI PLAN TO GET POINT ID
/// RETURNS PLAN TO GET ONE POINT
/// MAYBE HACK THIS TO GET OUR EXPLODE
void* prepare_plan(void)
{
  TupleDesc rdesc;  
  Oid *argtype;
  void *re, *pplan;
  char sql[75];
  //Link *p;
  //int patchid;
  /// make a connection via spi
  SPI_connect();
  ////elog(INFO, "Preparing plan");
  //if (res < 0)
  //{
  //  //elog(ERROR, "cannot connect to SPI");
  //} 
  /// execute the sql statement
  sprintf(sql, "SELECT id from multistar_g37_l8 where id = 0");
  SPI_execute(sql, true, 1);
  //if (res != SPI_OK_SELECT)
  //{
  //  //elog(ERROR, "select failed?");
  //}
  //if (SPI_processed != 1)
  //{
  //  //elog(ERROR, "there must be exactly 1 row");
  //}w
  rdesc = SPI_tuptable->tupdesc;
  ////elog(INFO, "Got redsc");
  argtype = (Oid *) palloc(sizeof(Oid));
  argtype[0] = SPI_gettypeid(rdesc, 1);
  argtype[1] = SPI_gettypeid(rdesc, 1);
  //elog(INFO,"PP Patch, point seems to be: %i, %i",argtype[0],argtype[1]);
  //pplan = SPI_prepare("SELECT * from (SELECT (gettin2(tree,points,stars)).* FROM multistar_g37_l8 WHERE id = $2) as f where f.id = $1 limit 1", 2, argtype);
  pplan = SPI_prepare("SELECT (getrow($1,tree,points,stars)).* FROM multistar_g37_l8 WHERE id = $2", 2, argtype); // Should be at least 5 times as fast on higher ids.
  re = SPI_saveplan(pplan); //SPI_keepplan?
///  SPI_finish();
  return re;
}

PG_FUNCTION_INFO_V1(testmulti);
Datum
testmulti(PG_FUNCTION_ARGS)
{
  void *pplan,*pplan2;
  Multistar **row;
  Vstar *s, *m;
  int32 point, patch;
  point = PG_GETARG_INT32(0);
  patch = PG_GETARG_INT32(1);
  pplan = prepare_plan();
  pplan2 = prepare_plan2();

  row[0] = fetch_row(pplan2,patch);
  //elog(INFO,"Got row");
  s = fetch_row_vstar(point,pplan,patch);
  //elog(INFO,"Got row");
  m = fetch_row_vstar2(point,row,pplan2,patch);
  //elog(INFO,"Got row");


 //elog(INFO,"XYZ %f %f %f %f %f %f",s->coords[0],m->coords[0],s->coords[1],m->coords[1],s->coords[2],m->coords[2]);
 //elog(INFO,"DIM %i %i",*(s->dim),*(m->dim));
 //elog(INFO,"STAR %i %i %i %i",s->star[0],m->star[0],s->star[1],m->star[1]);
  SPI_finish();
  PG_RETURN_INT32(1);
}

void* prepare_plan2(void)
/*  Function to get a row plan up
    when the current version doesn't work */
{
  TupleDesc rdesc;  
  Oid *argtype;
  void *re, *pplan;
  char sql[75];
  //Link *p;
  //int patchid;

  SPI_connect();

  sprintf(sql, "SELECT id from multistar_g37_l8 where id = 0");
  SPI_execute(sql, true, 1);
  rdesc = SPI_tuptable->tupdesc;
  ////elog(INFO, "Got redsc");
  argtype = (Oid *) palloc(sizeof(Oid));
  argtype[0] = SPI_gettypeid(rdesc, 1);
  //argtype[1] = SPI_gettypeid(rdesc, 1);
  //elog(INFO,"PP Patch, point seems to be: %i, %i",argtype[0],argtype[1]);
  //pplan = SPI_prepare("SELECT * from (SELECT (gettin2(tree,points,stars)).* FROM multistar_g37_l8 WHERE id = $2) as f where f.id = $1 limit 1", 2, argtype);
  pplan = SPI_prepare("SELECT id,tree,points,stars,bbox FROM multistar_g37_l8 WHERE id = $1", 1, argtype); // Should be at least 5 times as fast on higher ids.
  re = SPI_saveplan(pplan); //SPI_keepplan?
///  SPI_finish();
  return re;
}

Multistar* fetch_row(void* pplan2, int patchid)
/*  Fetches Multistar when needed */
{
  HeapTuple row;
  TupleDesc rdesc;
  ArrayType *tempa;
  bool isnull = false;
  Datum *parg;
  Multistar * multistar;
  SPITupleTable *tuptable;
  float4 *points;
  int i = 0;
  void *ptr;
  int size;
  parg = (Datum*) palloc(sizeof(Datum));
  multistar = (Multistar*) palloc(sizeof(Multistar));
  //nstars = 

  //elog(INFO,"FR Requested %i",patchid);

  parg[0] = Int32GetDatum(patchid);
  SPI_execute_plan(pplan2, parg, NULL, true, 1);
  row = SPI_tuptable->vals[0];
  rdesc = SPI_tuptable->tupdesc;

 //elog(INFO, "FR Actually got plan executed and it worked");
  multistar->patchid = patchid;
  //elog(INFO, "FR patchid"); 
  multistar->points = palloc(VARSIZE(DatumGetByteaP(SPI_getbinval(row, rdesc, 3, &isnull))));
  memcpy(multistar->points,((uint8_t*) DatumGetByteaP(SPI_getbinval(row, rdesc, 3, &isnull)))+4,VARSIZE(DatumGetByteaP(SPI_getbinval(row, rdesc, 3, &isnull))));
  //elog(INFO,"Points %f %f %f",multistar->points[0],multistar->points[1],multistar->points[2]);
  //multistar->points = (float4*) ( ( (uint8_t*) DatumGetByteaP(SPI_getbinval(row, rdesc, 3, &isnull)) ) + 4); // Y
  //elog(INFO, "FR points");

  multistar->nstars = (VARSIZE(DatumGetByteaP(SPI_getbinval(row, rdesc, 4, &isnull)))-VARHDRSZ)/sizeof(int32);
  //elog(INFO,"Nstars %i",multistar->nstars);
 //elog(INFO,"Multistars in total are %i",multistar->nstars);
/*
  size = VARSIZE(DatumGetByteaP(SPI_getbinval(row, rdesc, 4, &isnull)));//-VARHDRSZ;
  ptr = palloc(size);
  memcpy(ptr, DatumGetByteaP(SPI_getbinval(row, rdesc, 4, &isnull)), size);
*/
  multistar->stars = palloc(VARSIZE(DatumGetByteaP(SPI_getbinval(row, rdesc, 4, &isnull))));
  memcpy(multistar->stars,((uint8_t*) DatumGetByteaP(SPI_getbinval(row, rdesc, 4, &isnull)))+4,VARSIZE(DatumGetByteaP(SPI_getbinval(row, rdesc, 4, &isnull)))-4);
  //elog(INFO,"Stars %i %i %i",multistar->stars[0],multistar->stars[1],multistar->stars[2]);

  //multistar->stars = (int32*) ( ( (uint8_t*) DatumGetByteaP(SPI_getbinval(row, rdesc, 4, &isnull)) ) + 4); // Z

 //elog(INFO, "FR stars");
  tempa = DatumGetArrayTypeP(SPI_getbinval(row, rdesc, 2, &isnull)); // offsets
  //elog(INFO, "FR array");
  multistar->offsets = palloc(VARSIZE(tempa));
  memcpy(multistar->offsets,ARR_DATA_PTR(tempa),VARSIZE(tempa));
  //elog(INFO,"Offsets %i %i %i",multistar->offsets[0],multistar->offsets[1],multistar->offsets[2]);
  //multistar->offsets = (int*) ARR_DATA_PTR(tempa);
  
  //elog(INFO, "FR offsets %i",sizeof(BOX3D));
  multistar->bbox = palloc(sizeof(BOX3D));
  memcpy(multistar->bbox,DatumGetPointer(SPI_getbinval(row, rdesc, 5, &isnull)),sizeof(BOX3D));
  //multistar->bbox = (BOX3D*) DatumGetPointer(SPI_getbinval(row, rdesc, 5, &isnull));
  /*for (i=0;i<abs(multistar->offsets[0]);i++) {
    //elog(INFO,"%i",multistar->offsets[i]);
  }*/
  pfree(parg);
  tuptable = SPI_tuptable;
  SPI_freetuptable(tuptable);
  //elog(INFO,"FR Returning row %i",multistar->patchid);
  //elog(INFO,"FR Returning stars, %i %i %i",multistar->stars[0],multistar->stars[2],multistar->stars[3]);
  //elog(INFO,"Returning multistar %i",multistar->patchid);
  return multistar;
}


Vstar* fetch_row_vstar2(int id, Multistar** row, void* pplan_row, int patchid)
{
  Vstar *tempvstar;
  Link *p;
  int *starlen,len,l,i,a,b;
  int32 *results;
  int c = -1;

  starlen = palloc(sizeof(int32));
  //elog(INFO,"Got asked for patch %i point %i",patchid,id);


  if (id > 65536)
  {
    p = unstitch(id);
    id = p->point;
    patchid = p->patch;
    //elog(INFO,"FRV2 Unstitched to %i, %i",id,patchid);
  }

  /* Check if row exists in memory */
  for (i=0;i<3;i++) {
    if (patchid == row[i]->patchid) {
      c = i;
      break;
    }
  }

  /* Not optimally ordered 
  if (c == 1) {
    row[2] = row[0];
    row[0] = row[1];
    row[1] = row[2];
  }*/
  /* Else let's get it from disk */
  if (c == -1) {
    //pfree(row[2]->stars);
    //pfree(row[2]->points);
    //pfree(row[2]->bbox);
    //pfree(row[2]->offsets);
    //pfree(row[2]);
    row[2] = row[1];
    row[1] = row[0];
    //elog(INFO, "FRV2 Got wrong patchid req %i, is %i, new row here!",patchid,row[0]->patchid);
    row[0] = fetch_row(pplan_row,patchid);
    //memcpy(row[0],fetch_row(pplan_row,patchid),sizeof(Multistar));
    c = 0;
    //pfree(row[1]);
  }
  tempvstar = (Vstar*) palloc(sizeof(Vstar));
  len = abs(row[c]->offsets[0]);

  if (!(id > 0 && id <= len)) {
    elog(ERROR,"FRV2 Requested id not in bucket, strange?");
  }
  // END OF LIST
  else if (id == len) {
    //elog(INFO,"Seems we have a point at the end :O");
    b = VARSIZE(row[c]->stars)-VARHDRSZ;
    a = abs(row[c]->offsets[id]);
    //*starlen = (abs(b)-abs(a))/sizeof(int32);
    *starlen = row[c]->nstars-(abs(a)/sizeof(int32));
  }
  // USING OFFSET LIST
  else {
    a = abs(row[c]->offsets[id]);
    b = abs(row[c]->offsets[id+1]);
    *starlen = (abs(b)-abs(a))/sizeof(int32);
  }
  //elog(INFO,"Starlen %i",*starlen);

  tempvstar->id = id;
  tempvstar->flag = false;
  tempvstar->patchid = patchid;

  i = id -1;
  tempvstar->coords[0] = row[c]->points[3*i]; // X
  tempvstar->coords[1] = row[c]->points[3*i+1]; // Y
  tempvstar->coords[2] = row[c]->points[3*i+2]; // Z
  //elog(INFO,"FRV2 Returning xyz, %f %f %f",tempvstar->coords[0],tempvstar->coords[1],tempvstar->coords[2]);

  results = (int32*) palloc(sizeof(int32)*(*starlen));

  for (l=0;l<*starlen;l++) {
    results[l] = row[c]->stars[l+(a/sizeof(int32))];
    //elog(INFO, "FRV2 star[%i] is %i",l,row->stars[l+(a/sizeof(int32))]);
  }
  //elog(INFO,"FRV2 Returning star %i %i",results[0],results[1]);
  

  tempvstar->star = results;//row->stars[a/sizeof(int32)];
  tempvstar->dim = starlen;
  tempvstar->sid = stitch(patchid,id);
  return tempvstar;

}


/// FETCH A VSTAR FROM A ROW
/// USES PLAN TO GET TABLE ITSELF?
Vstar* fetch_row_vstar(int id, void* pplan, int patchid)
{
  HeapTuple row;
  TupleDesc rdesc;
  ArrayType *tempa;
  bool isnull = false;
  Datum *parg;
  Vstar *tempvstar;
  SPITupleTable *tuptable;
  Link *p;
  //int oldpatchid;


  parg = (Datum*) palloc(sizeof(Datum));
  tempvstar = (Vstar*) malloc(sizeof(Vstar));

  if (id > 65536)
  {
    p = unstitch(id);
    id = p->point;
    patchid = p->patch;
    //elog(INFO,"FRV Patch, point seems to be: %i, %i",id,patchid);
  }
  //elog(INFO,"FRV Requested %i",id);
  parg[0] = Int32GetDatum(id);
  parg[1] = Int32GetDatum(patchid);
  SPI_execute_plan(pplan, parg, NULL, true, 1);
  row = SPI_tuptable->vals[0];
  rdesc = SPI_tuptable->tupdesc;
  tempvstar->id = id;
  tempvstar->flag = false;
  tempvstar->patchid = patchid;
  tempvstar->sid = stitch(patchid,id);
  tempvstar->coords[0] = DatumGetFloat8(SPI_getbinval(row, rdesc, 2, &isnull)); // X
  tempvstar->coords[1] = DatumGetFloat8(SPI_getbinval(row, rdesc, 3, &isnull)); // Y
  tempvstar->coords[2] = DatumGetFloat8(SPI_getbinval(row, rdesc, 4, &isnull)); // Z
  tempa = DatumGetArrayTypeP(SPI_getbinval(row, rdesc, 5, &isnull)); // STAR
  tempvstar->star = (int*) ARR_DATA_PTR(tempa);
  tempvstar->dim = ARR_DIMS(tempa);
  //pfree(parg);
  //tuptable = SPI_tuptable;
  //SPI_freetuptable(tuptable);
  //elog(INFO,"FRV Returning vstar %i %i",tempvstar->id,tempvstar->patchid);
  return tempvstar;

}

/// RETURN POINT IN TRIANGLE
/// JUMP AND WALK ONLY, GOOD CANDIDATE
PG_FUNCTION_INFO_V1(point_location);
Datum point_location(PG_FUNCTION_ARGS)
{
  float8 dest[2];
  int startid;
  float8 power;
  Triangle *result;
  Vstar **tr;
  int patchid;
  int32 test;
  Link *p;
    uint8 *bytes;
  size_t bytes_size;
  bytea *wkb;
  size_t wkb_size;
  void *pplan2;

  //elog(INFO,"Got started with input.");


  dest[0] = PG_GETARG_FLOAT8(0); /// X
  dest[1] = PG_GETARG_FLOAT8(1); /// Y
  startid = PG_GETARG_INT32(2);  /// Start somewhere
  power   = PG_GETARG_FLOAT8(3);  /// ? power...
  patchid = PG_GETARG_INT32(4); // patchid
  SPI_connect();
  pplan2 = prepare_plan2();
  /////elog(WARNING, "spi connects");
  tr = jump_and_walk(dest, pplan2, startid, power, patchid);
  //elog(WARNING, "point_location returned tr %d %d %d", tr[0]->id, tr[1]->id, tr[2]->id);
  if (tr == NULL) {
    //elog(WARNING, "yes it's outside");
    SPI_finish();
    PG_RETURN_NULL();
  }
  else {
    SPI_finish();
    // COULD BE ACTUAL COORDINATES IN HERE
    bytes = trianglez_to_wkb(tr, &bytes_size);
    wkb_size = VARHDRSZ + bytes_size;
    wkb = palloc(wkb_size);
    memcpy(VARDATA(wkb), bytes, bytes_size);
    SET_VARSIZE(wkb, wkb_size);
    //elog(INFO,"Got wkb.");

    //pfree(bytes);

    SPI_finish();
    PG_RETURN_BYTEA_P(wkb);
  }
}

void* trianglez_to_wkb(Vstar **tr, size_t *wkbsize)
{
  /* With a given array of points and one triangle, we cast this triangle into a 
  well known binary form */

  uint32_t wkbtype = 1003; // POLYGONZ 1017; /* WKB TRIANGLEZ */
  uint32_t nump = 4; /* WKB TRIANGLEZ */
  uint32_t numr = 1; /* WKB TRIANGLEZ */

  size_t size = 1 + 4 + 4 + 4 + (4*8*3); /* endian + type + rings + num + 4xyz doubles */
  uint8_t *wkb, *ptr;
  Link *a, *b, *c;
  
  wkb = palloc(size);
  ptr = wkb;

  ptr[0] = machine_endian(); /* Endian flag */
  ptr += 1;

  memcpy(ptr, &wkbtype, 4); /* WKB type */
  ptr += 4;

  memcpy(ptr, &numr, 4); /* num of linear rings type */
  ptr += 4;

  memcpy(ptr, &nump, 4); /* num of points type */
  ptr += 4;

  /* triangle vertices are stitched as point, patch pair */

  memcpy(ptr, tr[0]->coords, 24); /* A */
  ptr += 24;

  memcpy(ptr, tr[1]->coords, 24); /* A */
  ptr += 24;

  memcpy(ptr, tr[2]->coords, 24); /* A */
  ptr += 24;

  memcpy(ptr, tr[0]->coords, 24); /* A */
  ptr += 24;
  
  if ( wkbsize ) *wkbsize = size;
  return wkb;
}

/// GET STARTING ID FROM buckets TABLE
/// SHOULD BE REPLACED OR DELETED
/// BY FUNCTION THAT JUST TAKES MIDDLE OF # POINTS
// SHOULD FIX FOR POINT OUTSIDE BBOX! Return -2
int get_starting_id(float8 dest[], int patchid)
{
  int res;
  char sql[400];
  double minx, maxx;
  int count;
  double min[2];
  double max[2];
  float len;
  float len2;
  float ratio;
  HeapTuple row;
  TupleDesc rdesc;
  bool isnull = false;

  SPI_connect();
  sprintf(sql, "SELECT ST_xMIN(ST_FORCE3DZ(bbox)) as minx, ST_xMAX(ST_FORCE3DZ(bbox)) as maxx,@tree[1] as count, ST_yMIN(ST_FORCE3DZ(bbox)) as miny, ST_yMAX(ST_FORCE3DZ(bbox)) as maxy FROM multistar_g37_l8 WHERE id = %i LIMIT 1",patchid);
  SPI_execute(sql, true, 1);
  if (SPI_processed != 1)
  {
    //elog(WARNING, "GSI Couldn't get patch.");
    return 1; // Patch always has 1 point
  }
  else
  {
    rdesc = SPI_tuptable->tupdesc;
    row = SPI_tuptable->vals[0];
    min[0] = DatumGetFloat8(SPI_getbinval(row, rdesc, 1, &isnull));
    max[0] = DatumGetFloat8(SPI_getbinval(row, rdesc, 2, &isnull));
    count = DatumGetInt32(SPI_getbinval(row, rdesc, 3, &isnull));
    min[1] = DatumGetFloat8(SPI_getbinval(row, rdesc, 4, &isnull));
    max[1] = DatumGetFloat8(SPI_getbinval(row, rdesc, 5, &isnull));
    //elog(INFO,"Found min,max,count,patchid %f %f %i %i",min[0],min[1],count,patchid);
  }

  if (! inside_rectangle(dest, min, max))
  { 
    //return -2; // SHOULD BE -2 BUT FOR LINE INTS, BETTER NOT
    return -2;
  }

  len = max[0]-min[0];
  len2 = abs(max[0]-dest[0]);
  if (len2 == 0) {
    // probably on intersection
    len2 = 0.1;
  }
  ratio = len2/len;
  res = ratio * count;
  //elog(INFO,"GSI Length,other length, ratio, res, count: %f %f %f %i %i",len,len2,ratio,res,count);
  //res = (1-((maxx-dest[0])/(maxx-minx)))*count; // X POSITION IN WIDTH TIMES COUNT
  ////elog(INFO,"Starting id found is %i with %f as input",res,dest[0]);
  if ((res <= 0) | (res >= count))
  {
    res = -1;
    //elog(INFO,"GSI returning -1");
  }
  //elog(INFO,"GSI Starting id found is %i",res);
  SPI_finish();
  return (int)res;

}

int get_starting_id_convex(float8 dest[], int patchid)
{
  int res;
  char sql[400];
  HeapTuple row;
  TupleDesc rdesc;
  bool isnull = false;
  char table[50];

  //SELECT * FROM (SELECT (gettin_convex(tree, points, stars)).* FROM multistar_g37_l8 WHERE id =0) as f WHERE x IS NOT NULL
  sprintf(table, "multistar_g37_l8");
  sprintf(sql, "SELECT id, @(x-%f)+ @(y-%f) as c FROM (SELECT (gettin2(tree,points,stars)).* FROM %s WHERE id = %i) as f WHERE star[1] = 0 ORDER BY c LIMIT 1",dest[0],dest[1],table,patchid);
  SPI_execute(sql, true, 1);
  if (SPI_processed != 1)
  {
    //elog(WARNING, "GSI Couldn't get patch.");
    return 0; // Patch always has 1 point
  }
  else
  {
    rdesc = SPI_tuptable->tupdesc;
    row = SPI_tuptable->vals[0];
    res = DatumGetInt32(SPI_getbinval(row, rdesc, 1, &isnull));
    //elog(INFO,"GSIC Found closest point on convex hull! %i",res);
  }

  return (int)res;
}


/// POINT LOCATION
/// THEN CALCULATES Z BASED ON WEIGHT
/// PROVIDED BY AREA OF TRIANGLES WHICH ARE SUBDIVISIONS
/// OF POINT IN TRIANGLE
PG_FUNCTION_INFO_V1(estimate_z);
Datum estimate_z(PG_FUNCTION_ARGS)
{
  float8 dest[2];
  int startid;
  Vstar **tr;
  //int res;
  double a0, a1, a2, at;
  double estim;
  int patchid;
  void *pplan2;

  dest[0] = PG_GETARG_FLOAT8(0);
  dest[1] = PG_GETARG_FLOAT8(1);
  startid = PG_GETARG_INT32(2); 
  patchid = PG_GETARG_INT32(3);
  
  SPI_connect();
  pplan2 = prepare_plan2();

  tr = jump_and_walk(dest, pplan2, startid, -1, patchid);
  if (tr == NULL) {
    SPI_finish();
    PG_RETURN_NULL();
  }
  else {
    
    a0 = area_triangle(dest, tr[1]->coords, tr[2]->coords);
    a1 = area_triangle(dest, tr[2]->coords, tr[0]->coords);
    a2 = area_triangle(dest, tr[0]->coords, tr[1]->coords);
    at = area_triangle(tr[0]->coords, tr[1]->coords, tr[2]->coords);
    estim = ( (a0 * tr[0]->coords[2]) + (a1 * tr[1]->coords[2]) + (a2 * tr[2]->coords[2]) ) / at; 
    //pfree(tr[0]);
    //pfree(tr[1]);
    //pfree(tr[2]);
    SPI_finish();
    PG_RETURN_FLOAT8(estim);
  }
}

/// POINT LOCATION AGAIN
/// THIS TIME WITHOUT STARTING INPUT OR POWER
PG_FUNCTION_INFO_V1(pl);
Datum pl(PG_FUNCTION_ARGS)
{
  float8 dest[2];
  Triangle *result;
  Vstar **tr,*n;
  int patchid;
  int id;
  int32 test,mode;
  Link *p;
  void *pplan2;
  uint8 *bytes;
  size_t bytes_size;
  bytea *wkb;
  size_t wkb_size;

  dest[0] = PG_GETARG_FLOAT8(0); // X
  dest[1] = PG_GETARG_FLOAT8(1); // Y
  patchid = PG_GETARG_INT32(2); // Patchid
  mode = PG_GETARG_INT32(3); // Mode
  SPI_connect();
  pplan2 = prepare_plan2();
  //elog(INFO,"PL Start");
  id = get_starting_id(dest, patchid); /// get starting point here

  //elog(INFO,"PL id %i",id);
  if (id == -2)
	  tr = NULL;
  else if (id == -1)
	  tr = jump_and_walk(dest, pplan2, 0, 0.2, patchid);
  else
	  tr = jump_and_walk(dest, pplan2, id, 0.2, patchid);
  if (tr == NULL) {
    SPI_finish();
    PG_RETURN_NULL();
  }
  else {
    elog(INFO,"Return bytea");
    bytes = trianglez_to_wkb(tr, &bytes_size);
    wkb_size = VARHDRSZ + bytes_size;
    wkb = palloc(wkb_size);
    memcpy(VARDATA(wkb), bytes, bytes_size);
    SET_VARSIZE(wkb, wkb_size);

    if (! mode == 0) {
      if (mode == 1){
        /* Aspect */
        Point *p[3];
        p[0] = palloc(sizeof(Point));
        p[1] = palloc(sizeof(Point));
        p[2] = palloc(sizeof(Point));

        Vec *n;
        float a;
        p[0]->x = (float4) tr[0]->coords[0];
        p[0]->y = (float4) tr[0]->coords[1];
        p[0]->z = (float4) tr[0]->coords[2];
        p[1]->x = (float4) tr[1]->coords[0];
        p[1]->y = (float4) tr[1]->coords[1];
        p[1]->z = (float4) tr[1]->coords[2];
        p[2]->x = (float4) tr[2]->coords[0];
        p[2]->y = (float4) tr[2]->coords[1];
        p[2]->z = (float4) tr[2]->coords[2];
        n = normal(p);
        a = aspect(n);
        elog(INFO,"Aspect: %f",a);  
      }
      else if (mode == 2){
        /* Slope */
        Point *p[3];
        p[0] = palloc(sizeof(Point));
        p[1] = palloc(sizeof(Point));
        p[2] = palloc(sizeof(Point));

        Vec *n;
        float s;
        p[0]->x = (float4) tr[0]->coords[0];
        p[0]->y = (float4) tr[0]->coords[1];
        p[0]->z = (float4) tr[0]->coords[2];
        p[1]->x = (float4) tr[1]->coords[0];
        p[1]->y = (float4) tr[1]->coords[1];
        p[1]->z = (float4) tr[1]->coords[2];
        p[2]->x = (float4) tr[2]->coords[0];
        p[2]->y = (float4) tr[2]->coords[1];
        p[2]->z = (float4) tr[2]->coords[2];
        n = normal(p);
        s = slope(n);
        elog(INFO,"Slope: %f",s);  
      }
      else if (mode == 3){
        /* Local min of first point */
        Multistar **row;
        int i;  
        bool s = true;
        row = (Multistar**) palloc(sizeof(Multistar*)); // Complete fucking magic
        row[0] = (Multistar*) palloc(sizeof(Multistar)); // Actually storing list of 
        row[1] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
        row[2] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
        for (i=0;i<*tr[0]->dim;i++) {
          n = fetch_row_vstar2(tr[0]->star[i],row,pplan2,patchid);
          if (n->coords[2] < tr[0]->coords[2]) {
            elog(INFO,"False");
            s = false;
            break;
          }
        }
        if (s) { elog(INFO,"True"); }
      }
      else if (mode == 4){
        /* Local max of first point */
        Multistar **row;
        int i;
        bool s = true;
        row = (Multistar**) palloc(sizeof(Multistar*)); // Complete fucking magic
        row[0] = (Multistar*) palloc(sizeof(Multistar)); // Actually storing list of 
        row[1] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
        row[2] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
        for (i=0;i<*tr[0]->dim;i++) {
          n = fetch_row_vstar2(tr[0]->star[i],row,pplan2,patchid);
          if (n->coords[2] > tr[0]->coords[2]) {
            elog(INFO,"False");
            s = false;
            break;
          }
        }
        if (s) {elog(INFO,"True");}
      }
      else if (mode == 5){
        /* Degree of first point */
        elog(INFO,"Degree: %i",*tr[0]->dim);
      }
    }
    //pfree(bytes);
    //SPI_finish();
    PG_RETURN_BYTEA_P(wkb);
  }
}

/*
Vstar* walk_convex(int start, float8 a, float8 b, int patchid)
// Actually trying to walk around the border, intersecting line segments
// If we find matching line segment, rejoice. Use triangle as return for 
// walking straight.
{
    return 1;
}
*/

/// SAME AS PREVIOUS
/// NOW RETURNS TOTAL INTERSECTS
PG_FUNCTION_INFO_V1(profile_count_intersections);
Datum profile_count_intersections(PG_FUNCTION_ARGS)
{
  float8 a[2], b[2];
  Vstar **tr;
  //int res;
  int total_intersects;
  void *pplan;
  int bucketid;
  int patchid;
  bytea *wkb;
  void *bytes;
  size_t size_o;
  size_t bytes_size;// = (size_t*) palloc(sizeof(size_t));
  size_t wkb_size;

  a[0] = PG_GETARG_FLOAT8(0);
  a[1] = PG_GETARG_FLOAT8(1);
  b[0] = PG_GETARG_FLOAT8(2);
  b[1] = PG_GETARG_FLOAT8(3);
  patchid = PG_GETARG_INT32(4);

  SPI_connect();
  pplan = prepare_plan2();

  bucketid = get_starting_id(a, patchid);
  elog(INFO,"PCI Starting id: %i",bucketid);
  if (bucketid == -2)
  {
    // POINT OUTSIDE BBOX.
    // LETS USE CONVEX HULL! 
    // YOU SHOULD ACTUALLY USE extended gettin function!
    // WHICH FILTERS ON THE FLY, NOT AFTERWARDS (by -offset)
    tr = get_starting_triangle_convex(a,b,patchid);
    //elog(INFO, "PCI Found a triange by convex.");
  }

  else {
    tr = jump_and_walk(a, pplan, bucketid, -1, patchid);
   //elog(INFO, "PCI Found a triangle by jumping and walking.");
    if (tr == NULL) {
      tr = get_starting_triangle_convex(a,b,patchid);
     //elog(INFO, "PCI Found a triange by convex after JW found nothing.");
    }
  }
  //elog(INFO,"PCI Starting id : %i",bucketid);
  if (tr == NULL) {
    //SPI_finish();
   //elog(INFO,"No starting triangle, returning NULL");
    PG_RETURN_NULL();
  }
  //pplan = prepare_plan();  // single
  //SPI_connect();
  //pplan = prepare_plan2();
  //elog(INFO,"Actually got plan!");

  tr[0]->flag = false;
  tr[1]->flag = false;
  tr[2]->flag = false;
  elog(INFO,"Triangle ids %i %i %i",tr[0]->id,tr[1]->id,tr[2]->id);

  //total_intersects = walk_straight_and_report_intersections(tr, a, b, pplan, false, patchid);
  bytes = walk_straight_and_report_intersections(tr, a, b, pplan, false, patchid, &bytes_size);

  /* How to make a LineStringZ

  <linestringz binary representation> ::= 
      <byte order> <wkblinestringz> <num> <wkbpointz binary>

  
  */
  elog(INFO,"Got %i",bytes_size);
  //SPI_finish();

  wkb_size = VARHDRSZ + bytes_size;
  wkb = palloc(wkb_size);
  memcpy(VARDATA(wkb), bytes, bytes_size);
  SET_VARSIZE(wkb, wkb_size);
  //elog(INFO,"Got wkb.");

  pfree(bytes);
  //elog(INFO,"Returning.");




  elog(INFO,"Hi");
  SPI_freeplan(pplan);
  //elog(INFO,"PCI Completely done \n");
  ////elog(INFO,"====================\n");
  //PG_RETURN_INT32(total_intersects);
  elog(INFO,"Hi2");
  PG_RETURN_BYTEA_P(wkb);
}

void* walk_straight_and_report_intersections(Vstar **tr, float8 a[], float8 b[], void *pplan, bool justcount, int patchid, size_t *wkbsize)
{
  bool outside = false;
  Vstar *tempvstar;
  Link *p;
  int prev;
  double x[3];
  int total_intersects = 0;
  Multistar **multistar;// = palloc(sizeof(Multistar));
  
  ///////////////////////////////////////////

 //elog(INFO,"WSRI Got here");
  uint32_t wkbtype = 1002;
  size_t size = 1 + 4 + 4 +(24 * 100000); /* endian + type + numt + numt*size of TriangleZ */
  bytea *wkb;
  uint8_t *wkbs, *ptr;
  uint32_t wkb_size;
  wkbs = palloc(size);
  ptr = wkbs;

  ptr[0] = machine_endian(); /* Endian flag */
  ptr += 1;

  memcpy(ptr, &wkbtype, 4); /* WKB type */
  ptr += 4;
  ptr += 4; // npoints

  ///////////////////////////////////////////
  //SPI_connect();
 //elog(INFO,"WSRI Past preparation");

  //parg = (Datum*) palloc(sizeof(Datum));
  multistar = (Multistar**) palloc(sizeof(Multistar*)*3); // Complete fucking magic
  multistar[0] = (Multistar*) palloc(sizeof(Multistar)); // Actually storing list of 
  multistar[1] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
  multistar[2] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
  multistar[0] = fetch_row(pplan,patchid);

  //elog(INFO,"WSRI Allocated multistars");
  //elog(WARNING, "WSR \n-----INTERSECTIONS-----\n");

  if (intersect(tr[2]->coords, tr[0]->coords, a, b) == true) {
    //elog(INFO,"WSR Rotating");
    tempvstar = tr[0];
    tr[0] = tr[1];
    tr[1] = tr[2];
    tr[2] = tempvstar;
  }
 //elog(INFO,"WSR Starting to loop.");
  while (true) {
    if (intersect(tr[0]->coords, tr[1]->coords, a, b) == false) { /* ==0 OR == 1 */
      if (intersect(tr[1]->coords, tr[2]->coords, a, b) == false) { /* ==0 OR == 1 */
        /* found! */
        //elog(INFO,"FOUND! %i %i %i",tr[0]->id, tr[1]->id, tr[2]->id);
        break;
      }
      else {
        /* jump: (1, prev, 2) */
        if (tr[1]->patchid != tr[2]->patchid) {
          prev = find_previous_in_star(tr[1]->star, *(tr[1]->dim), stitch(tr[2]->patchid,tr[2]->id));
          patchid = tr[1]->patchid;
        }
        else {
          prev = find_previous_in_star(tr[1]->star, *(tr[1]->dim), tr[2]->id);
          patchid = tr[1]->patchid;
        }
        if (prev == 0) {
          outside = true;
          //elog(INFO,"Outside! %i %i %i",tr[0]->id, tr[1]->id, tr[2]->id);
          break;
        }/*
        if (prev > 65536)
        {
          p = unstitch(prev);
          prev = p->point;
          patchid = p->patch;
        }*/
        tempvstar = tr[0];
        tr[0] = tr[1];
        //tr[1] = tempvstar;
        //tr[1]->id = prev;
        //tr[1]->flag = false;
        //patchid = tr[1]->patchid;
        //tr[1] = fetch_row_vstar(prev, pplan,patchid);
        tr[1] = fetch_row_vstar2(prev, multistar, pplan,patchid);
        patchid = tr[1]->patchid;
        //tr[1]->flag = false;
      }
    }
    else {
      /* jump: (0, prev, 1) */
      if (tr[0]->patchid != tr[1]->patchid) {
        prev = find_previous_in_star(tr[0]->star, *(tr[0]->dim), stitch(tr[1]->patchid,tr[1]->id));
        patchid = tr[0]->patchid;
      }
      else {
        prev = find_previous_in_star(tr[0]->star, *(tr[0]->dim), tr[1]->id);
        patchid = tr[0]->patchid;
      }
      if (prev == 0) {
        outside = true;
        //elog(INFO,"Outside! %i %i %i",tr[0]->id, tr[1]->id, tr[2]->id);
        break;
      }/*
      if (prev > 65536)
      {
        p = unstitch(prev);
        prev = p->point;
        patchid = p->patch;
      }*/
      tempvstar = tr[1];
      tr[1] = tr[2];
      tr[2] = tempvstar;
      //tr[1]->id = prev;
      //tr[1]->flag = false;
      //tr[1] = fetch_row_vstar(prev, pplan,patchid);
      tr[1] = fetch_row_vstar2(prev, multistar, pplan,patchid);
      patchid = tr[1]->patchid; // previously disabled?
      //tr[1]->flag = false;
    }
    /* report here the intersection */
    //elog(INFO, "LINESTRING(%f %f, %f %f)",tr[0]->coords[0],tr[0]->coords[1], tr[2]->coords[0],tr[2]->coords[1]);
    /// calculate the z of intersection
    /*if (justcount == false)
    {*/
      get_intersection_z(a, b, tr[0]->coords, tr[2]->coords, x);
      memcpy(ptr,x,24);
      ptr += 24;

      

      //elog(WARNING, "intersect edge (%d--%d), and x is %f %f %f", tr[0]->id, tr[2]->id, x[0], x[1], x[2]);
      
      
      /* Gather all x here and copy them to object every time. Use /*

    }*/
    total_intersects++;
    //elog(WARNING, "%f %f %f", x[0], x[1], x[2]);
  }
  //elog(WARNING, "WSRI Total intersections is %d", total_intersects);
  
  memcpy(wkbs+5, &total_intersects, 4); /* num of points */

  size = 1 + 4 + 4 + ( total_intersects * 24 );
  if ( wkbsize ) *wkbsize = size;
  //SPI_finish();
  //wkb_size = VARHDRSZ + size;
  //wkb = palloc(wkb_size);
  //memcpy(VARDATA(wkb), wkbs, size);
  //SET_VARSIZE(wkb, wkb_size);
  elog(INFO,"Got wkb.");

  pfree(multistar[0]);
  pfree(multistar[1]);
  pfree(multistar[2]);
  pfree(multistar);
  //SPI_freeplan(pplan);
  elog(INFO,"Returning with %i ints, total size %i",total_intersects,size);
  return wkbs;
  //return total_intersects;
}


/* *************************************** */
void walk_straight(Vstar **tr, Multistar **multistar, float8 a[], float8 b[], float8 bboxmin[], float8 bboxmax[], void *pplan, int patchid)
{
  Datum *parg;
  SPITupleTable *tuptable;
  Vstar *tempvstar;
  Vstar *vstar1;
  int prev;
  bool outside;
  //Multistar **multistar;

  multistar = (Multistar**) palloc(sizeof(Multistar*)); // Complete fucking magic
  multistar[0] = (Multistar*) palloc(sizeof(Multistar)); // Actually storing list of 
  multistar[1] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
  multistar[2] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.

  elog(INFO, "WS Test");
  multistar[0] = fetch_row(pplan,patchid);
  
  elog(INFO,"WS Got rows");
  
  if (intersect(tr[2]->coords, tr[0]->coords, a, b) == true) {
    tempvstar = tr[0];
    tr[0] = tr[1];
    tr[1] = tr[2];
    tr[2] = tempvstar;
  }

  while (true) {
    if (intersect(tr[0]->coords, tr[1]->coords, a, b) == false) { /* ==0 OR == 1 */
      if (intersect(tr[1]->coords, tr[2]->coords, a, b) == false) { /* ==0 OR == 1 */
        /* found! */
        break;
      }
      else {
        /* jump: (1, prev, 2) */ 
        //prev = find_previous_in_star(tr[1]->star, *(tr[1]->dim), tr[2]->id);
        if (tr[0]->patchid != tr[1]->patchid) {
          prev = find_previous_in_star(tr[1]->star, *(tr[1]->dim), stitch(tr[2]->patchid,tr[2]->id));
          patchid = tr[1]->patchid;
        }
        else {
          prev = find_previous_in_star(tr[1]->star, *(tr[1]->dim), tr[2]->id);
          patchid = tr[1]->patchid;
        }
        if (prev == 0) {
          outside = true;
          break;
        }
        tempvstar = tr[0];
        tr[0] = tr[1];
        tr[1] = tempvstar;
        tr[1] = fetch_row_vstar2(prev, multistar, pplan, patchid);

        if (inside_rectangle(tr[2]->coords, bboxmin, bboxmax) == true) {
          HASH_FIND_INT(htable, &(tr[2]->sid), vstar1);
          if (vstar1 == NULL) {
            vstar1 = (Vstar*) malloc(sizeof(Vstar)); 
            *vstar1 = *(tr[2]);
            HASH_ADD_INT(htable, sid, vstar1);
          }
        } 
        tuptable = SPI_tuptable;
        SPI_freetuptable(tuptable);
      }
    }
    else {
      /* jump: (0, prev, 1) */
      //prev = find_previous_in_star(tr[0]->star, *(tr[0]->dim), tr[1]->id);
      if (tr[0]->patchid != tr[1]->patchid) {
        prev = find_previous_in_star(tr[0]->star, *(tr[0]->dim), stitch(tr[1]->patchid,tr[1]->id));
        patchid = tr[0]->patchid;
      }
      else {
        prev = find_previous_in_star(tr[0]->star, *(tr[0]->dim), tr[1]->id);
        patchid = tr[0]->patchid;
      }
      if (prev == 0) {
        outside = true;
        break;
      }
      tempvstar = tr[1];
      tr[1] = tr[2];
      tr[2] = tempvstar;
      tr[1] = fetch_row_vstar2(prev, multistar, pplan, patchid);
      if (inside_rectangle(tr[2]->coords, bboxmin, bboxmax) == true) {
        HASH_FIND_INT(htable, &(tr[2]->sid), vstar1);
        if (vstar1 == NULL) {
          vstar1 = (Vstar*) malloc(sizeof(Vstar)); 
          *vstar1 = *(tr[2]);
          HASH_ADD_INT(htable, sid, vstar1);
        }
      } 
      tuptable = SPI_tuptable;
      SPI_freetuptable(tuptable);
      //pfree(multistar[2]);
      //pfree(multistar[1]);
      //pfree(multistar[0]);
      //pfree(multistar);
    }
  } 
  elog(INFO,"WS Done");
}


/* *************************************** */
Vstar** get_starting_triangle_convex(float8 a[], float8 b[], int patchid)
{
  HeapTuple row;
  TupleDesc rdesc;
  bool isnull = false;
  Vstar **tr;
  Vstar *tempvstar;
  void *pplan;
  int pointid,i,o,or,ii=0;
  char sql[200];

  pplan = prepare_plan(); 

  //SELECT * FROM (SELECT (gettin_convex(tree, points, stars)).* FROM multistar_g37_l8 WHERE id =0) as f WHERE x IS NOT NULL
  sprintf(sql, "SELECT id, @(x-%f)+ @(y-%f) as c FROM (SELECT (gettin2(tree,points,stars)).* FROM multistar_g37_l8 WHERE id = %i) as f WHERE star[1] = 0 ORDER BY c LIMIT 1",a[0],a[1],patchid);
  SPI_execute(sql, true, 1);
  if (SPI_processed != 1)
  {
    //elog(WARNING, "GSIC Couldn't get patch.");
    return NULL; // Patch always has 1 point
  }
  else
  {
    rdesc = SPI_tuptable->tupdesc;
    row = SPI_tuptable->vals[0];
    pointid = DatumGetInt32(SPI_getbinval(row, rdesc, 1, &isnull));
    //elog(INFO,"GSIC Found closest point on convex hull! %i",pointid);
  }
  /* Now we've got an id of a nice point, hopefully somewhat close to our starting triangle.*/
  tr = (Vstar**) palloc(sizeof(Vstar*)); // Complete fucking magic
  tr[0] = (Vstar*) palloc(sizeof(Vstar)); // Actually storing list of 
  tr[1] = (Vstar*) palloc(sizeof(Vstar)); // mem locations to stars.
  tr[2] = (Vstar*) palloc(sizeof(Vstar)); // mem locations to stars.
  tempvstar = fetch_row_vstar(pointid,pplan,patchid);
  or = orient2d(a,b,tempvstar->coords);
  //elog(INFO,"GSIC Orientation is %i and point is %i",or,tempvstar->id);
  if (or == 1) {
    tr[0] = tempvstar;
    while (true) {
      i = tr[0]->star[1]; //# 0 is first, than we've got counter clock wise, coming from the left
      tr[2] = tr[0]; // move one back
      tr[0] = fetch_row_vstar(i,pplan,patchid);
      o = orient2d(a,b,tr[0]->coords);
      //elog(INFO,"GSIC is %i and it orients %i",i,o);
      if (o != or) {
        i = find_previous_in_star(tr[0]->star,*(tr[0]->dim),tr[2]->id);
        tr[1] = fetch_row_vstar(i,pplan,patchid);
        //elog(INFO,"GSIC Found starting triangle! %i %i %i",tr[0]->id, tr[1]->id,tr[2]->id);
        break;
      }
      if (i == pointid || ii > 500) {
        return NULL;
      }
      ii++;
    }
  }
  else if (or == -1) {
    tr[2] = tempvstar;
    while (true) {
      i = find_previous_in_star(tr[2]->star,*(tr[2]->dim),0); //# 0 is first, than we go clockwise, coming from the right
      tr[0] = tr[2];
      tr[2] = fetch_row_vstar(i,pplan,patchid);
      o = orient2d(a,b,tr[2]->coords);
      //elog(INFO,"GSIC is %i and it orients %i",i,o);
      if (o != or) {
        i = find_previous_in_star(tr[0]->star,*(tr[0]->dim),tr[2]->id);
        tr[1] = fetch_row_vstar(i,pplan,patchid);
        //elog(INFO,"GSIC Found starting triangle! %i %i %i",tr[0]->id, tr[1]->id,tr[2]->id);
        break;
      }
      if (i == pointid || ii > 500) {
        return NULL;
      }
      ii++;
    }
  }
  else {
    //elog(INFO,"GSIC Found nothing!");
    return NULL;
  }
  return tr;
}


Vstar** jump_and_walk(float8 dest[], void* pplan2, int startid, float8 power, int patchid)
{
  float8 curpt[2];
  HeapTuple row;
  TupleDesc rdesc;
  bool isnull = false;
  Vstar **tr;
  Vstar *tempvstar;
  int i, prev;
  int totalvisited = 0;
  void *pplan;
  Datum *parg;
  int dist;
  SPITupleTable *tuptable;
  bool outside = false;
  int nopts, samples, tempid;
  double mindist;
  Multistar **multistar;
  int a = 0; // switched 1 ++
  int b = 1; // switched 2 ++
  int c = 1; // switched 0 --
  int d = 2; // switched 1 --

  elog(INFO,"JW Startid %i patchid %i",startid,patchid);

  //SPI_connect();
  pplan = prepare_plan(); 
  //pplan2 = prepare_plan2();
  parg = (Datum*) palloc(sizeof(Datum));

  elog(INFO,"JW Startid %i patchid %i",startid,patchid);
  
  multistar = (Multistar**) palloc(sizeof(Multistar*)); // Complete fucking magic
  multistar[0] = (Multistar*) palloc(sizeof(Multistar)); // Actually storing list of 
  multistar[1] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
  multistar[2] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
  multistar[0] = fetch_row(pplan2,patchid);

  if (startid == 0) {
    //elog(WARNING, "JW #of points %i", no_points(patchid));
    //elog(WARNING, "%d %d", RAND_MAX, rand() % no_points(patchid));
    nopts = no_points(patchid);
    if (power < 0) 
      samples = pow(nopts, 0.25);
    else
      samples = pow(nopts, power);
    ////elog(WARNING, "JW #of samples checked: %d", samples);
    mindist = 1e8;
    /* finding a good starting vertex */
    for (i = 0; i < samples; i++) {
      parg[0] = Int32GetDatum(rand() % nopts);
      parg[1] = Int32GetDatum(patchid);
      SPI_execute_plan(pplan, parg, NULL, true, 1);
      row = SPI_tuptable->vals[0];
      rdesc = SPI_tuptable->tupdesc;
      curpt[0] = DatumGetFloat8(SPI_getbinval(row, rdesc, 2, &isnull));
      curpt[1] = DatumGetFloat8(SPI_getbinval(row, rdesc, 3, &isnull));
      if (distance_sqr(curpt, dest) < mindist) {
        mindist = distance_sqr(curpt, dest);
        startid =  DatumGetInt32(SPI_getbinval(row, rdesc, 1, &isnull));
      }
      tuptable = SPI_tuptable;
      SPI_freetuptable(tuptable); 
    }
    //elog(WARNING, "JW closest pt after random sampling is: %d", startid);
    /* getting a bit closer... */
    for (i = 0; i < 5; i++) {
      tempid = startid - (nopts * 0.1) + (rand() % (int)(nopts * 0.2));
      ////elog(WARNING, "--testing id: %d",tempid);

      parg[0] = Int32GetDatum(tempid);
      parg[1] = Int32GetDatum(patchid);
      SPI_execute_plan(pplan, parg, NULL, true, 1);
      if (SPI_processed != 0) {
        row = SPI_tuptable->vals[0];
        rdesc = SPI_tuptable->tupdesc;
        curpt[0] = DatumGetFloat8(SPI_getbinval(row, rdesc, 2, &isnull));
        curpt[1] = DatumGetFloat8(SPI_getbinval(row, rdesc, 3, &isnull));
        if (distance_sqr(curpt, dest) < mindist) {
          mindist = distance_sqr(curpt, dest);
          startid =  DatumGetInt32(SPI_getbinval(row, rdesc, 1, &isnull));
          ////elog(WARNING, "find one closer!");
        }
        tuptable = SPI_tuptable;
        SPI_freetuptable(tuptable);
      }
    }
  }
  //elog(WARNING, "JW closest pt to start is %d", startid);

  //parg[0] = Int32GetDatum(startid);
  //parg[1] = Int32GetDatum(patchid);
  //SPI_execute_plan(pplan, parg, NULL, true, 1);
  ////elog(INFO, "JW Got plan executed");

  row = SPI_tuptable->vals[0];
  rdesc = SPI_tuptable->tupdesc;

  tr = (Vstar**) palloc(sizeof(Vstar*)); // Complete fucking magic
  tr[0] = (Vstar*) palloc(sizeof(Vstar)); // Actually storing list of 
  tr[1] = (Vstar*) palloc(sizeof(Vstar)); // mem locations to stars.
  tr[2] = NULL;
  elog(INFO,"Alloced Vstar");
  ///-- tr[0]
  tr[0] = fetch_row_vstar2(startid, multistar, pplan2, patchid);
  elog(INFO,"Triangle %i",tr[0]->sid);

  /* check the distance */
  dist = sqrt(distance_sqr(dest, tr[0]->coords));
  //elog(WARNING, "JW start distance is %f", dist);

  /*  find the first triangle to start the walk */
  /*  Looping over star, do not switch context by changing 
      patchid. */
  for (i = 0; i < *(tr[0]->dim); i++) {
    if (tr[0]->star[i] != 0) {
      //parg[0] = Int32GetDatum(tr[0]->star[i]);
      //parg[1] = Int32GetDatum(patchid);
      //elog(INFO, "JW Trying vertex %i",tr[0]->star[i]);
      //SPI_execute_plan(pplan, parg, NULL, true, 1);
      //elog(WARNING, "JW id star %d", tr[0]->star[i]);
      /*row = SPI_tuptable->vals[0];
      rdesc = SPI_tuptable->tupdesc;
      tr[1]->id = DatumGetInt32(SPI_getbinval(row, rdesc, 1, &isnull));
      tr[1]->coords[0] = DatumGetFloat8(SPI_getbinval(row, rdesc, 2, &isnull));
      tr[1]->coords[1] = DatumGetFloat8(SPI_getbinval(row, rdesc, 3, &isnull));
      tr[1]->coords[2] = DatumGetFloat8(SPI_getbinval(row, rdesc, 4, &isnull));
      tempa = DatumGetArrayTypeP(SPI_getbinval(row, rdesc, 5, &isnull));
      tr[1]->dim = ARR_DIMS(tempa);
      tr[1]->star = (int*) ARR_DATA_PTR(tempa);
      tr[1]->patchid = patchid;*/
      tr[1] = fetch_row_vstar2(tr[0]->star[i],multistar,pplan2,patchid);
      totalvisited++;

      //tuptable = SPI_tuptable;
      //SPI_freetuptable(tuptable);
      //elog(INFO, "Looping over star of first point, finding edge that fits.");

      if (orient2d(tr[0]->coords, tr[1]->coords, dest) == -1) { 
        //elog(INFO, "Orient true");

        /// find previous vertex in star(startid)
        prev = find_previous_in_star(tr[0]->star, *(tr[0]->dim), tr[1]->id);

        if (tr[0]->patchid != tr[1]->patchid) {
          prev = find_previous_in_star(tr[0]->star, *(tr[0]->dim), stitch(tr[1]->patchid,tr[1]->id));
          patchid = tr[1]->patchid;
        }
        else {
          prev = find_previous_in_star(tr[0]->star, *(tr[0]->dim), tr[1]->id);
          patchid = tr[1]->patchid;
        }

        if (prev == 0) { // We've got a that probably outside, so selecting next best thing.
         //elog(INFO,"Prev is 0");
          prev = find_next_in_star(tr[0]->star, *(tr[0]->dim), tr[1]->id);
        }
        /* vstar 2 */
        //elog(INFO,"Previous is %i, with %i as start and %i as correct.",prev,tr[0]->id,tr[1]->id);
        tr[2] = tr[1];
        /* vstar 1 */
        tr[1] = (Vstar*) palloc(sizeof(Vstar)); 
        tr[1] = fetch_row_vstar2(prev, multistar, pplan2, patchid);

        //tuptable = SPI_tuptable;
        //SPI_freetuptable(tuptable);
        //elog(INFO, "JW Found something");
        break;
      }
    }
  }
  //if (!(tr[0]->id != tr[1]->id && tr[1]->id != tr[2]->id && tr[2]->id != tr[0]->id)) {
    //elog(WARNING,"HORRIBLE THINGS WILL HAPPEN WITH THIS NON TRIANGLE");
  //  return NULL;
  //}

  /*//elog(WARNING, "JW the first triangle is (%d, %d, %d)", tr[0]->id, tr[1]->id, tr[2]->id);
  //elog(WARNING, "JW tr visited %d", totalvisited);
  //elog(WARNING, "JW --- (%f, %f, %f)", tr[0]->coords[0], tr[1]->coords[0], tr[2]->coords[0]);
  //elog(WARNING, "JW --- (%f, %f, %f)", tr[0]->coords[1], tr[1]->coords[1], tr[2]->coords[1]);*/

  /*  From here on we could walk into another patch... So we should switch patchid for that...
      As soon as we find a new point, we should/could wait? 
  */
  elog(INFO,"Started looping %i %i %i",tr[0]->id,tr[1]->id,tr[2]->id);


  while (true) {
    //elog(INFO,"JW looping with Triangle of %i %i %i",tr[0]->id,tr[1]->id,tr[2]->id);
    /* The following does two orientation tests, which should be done in a random order
    if we're walking a constrained Delaunay dataset. 
    */

    
    /*
    if (a == 0){ a = 1; } else {a = 0;}
    if (b == 1){ a = 2; } else {a = 1;}
    if (c == 1){ a = 0; } else {a = 1;}
    if (d == 2){ a = 1; } else {a = 2;}
    */

    if (orient2d(tr[a]->coords, tr[b]->coords, dest) != -1) { /* ==0 OR == 1 */
      //elog(INFO,"0");

      if (orient2d(tr[c]->coords, tr[d]->coords, dest) != -1) { /* ==0 OR == 1 */
        /* found! */
        //elog(INFO, "JW Found!");
        break;
      }
      else {
        /* jump: (1, prev, 2) */ 
        if (tr[c]->patchid != tr[d]->patchid) {
          prev = find_previous_in_star(tr[c]->star, *(tr[c]->dim), stitch(tr[d]->patchid,tr[d]->id));
          patchid = tr[c]->patchid;
        }
        else {
          prev = find_previous_in_star(tr[c]->star, *(tr[c]->dim), tr[d]->id);
          patchid = tr[c]->patchid;
        }
        if (prev == 0) {
          outside = true;
          break;
        }
        tempvstar = tr[a];
        tr[a] = tr[b];
        tr[b] = tempvstar;
        tr[b] = fetch_row_vstar2(prev, multistar, pplan2, patchid);
        //elog(INFO,"Fetched prev %i in patch %i",prev,patchid);
        //elog(WARNING, "tr %i %i %i", tr[0]->sid, tr[1]->sid, tr[2]->sid);

        //tuptable = SPI_tuptable;
        //SPI_freetuptable(tuptable);
        //elog(WARNING, "tr %d %d %d", tr[0]->id, tr[1]->id, tr[2]->id);
      }
    }
    else {
      /* jump: (0, prev, 1) */
      if (tr[a]->patchid != tr[b]->patchid) {
        prev = find_previous_in_star(tr[a]->star, *(tr[a]->dim), stitch(tr[b]->patchid,tr[b]->id));
        patchid = tr[a]->patchid;
      }
      else {
        prev = find_previous_in_star(tr[a]->star, *(tr[a]->dim), tr[b]->id);
        patchid = tr[a]->patchid;
      }
      if (prev == 0) {
        outside = true;
        break;
      }
      tempvstar = tr[c];
      tr[c] = tr[d];
      tr[d] = tempvstar;
      tr[c] = fetch_row_vstar2(prev, multistar, pplan2, patchid);
      //elog(INFO,"Fetched prev %i in patch %i",prev,patchid);
      //tuptable = SPI_tuptable;
      //SPI_freetuptable(tuptable);
      //elog(WARNING, "tr %i %i %i", tr[0]->sid, tr[1]->sid, tr[2]->sid);
    }
    //totalvisited++;
    //if (totalvisited % 10 == 0)
    //elog(WARNING, "tr visited %d", totalvisited);
  } 
  if (outside == true) {
    //elog(WARNING, "JW outside");
    pfree(tr[0]);
    pfree(tr[1]);
    pfree(tr[2]);
    SPI_freeplan(pplan);
    //SPI_freeplan(pplan2);
    //SPI_finish();
    return NULL;
  }
  else {
    ////elog(WARNING, "the found triangle is (%d, %d, %d)", tr[0]->id, tr[1]->id, tr[2]->id);
    elog(WARNING, "# of triangles visited is %d", totalvisited);
    SPI_freeplan(pplan);
    //SPI_freeplan(pplan2);
    //SPI_finish();
    //elog(INFO,"JW Found Triangle of %i %i %i",tr[0]->id,tr[1]->id,tr[2]->id);
    return tr;
  }
}
/*    result = (Triangle*) palloc(sizeof(Triangle));
    result->v0 = tr[0]->id;
    result->v1 = tr[1]->id;
    result->v2 = tr[2]->id;
    SPI_freeplan(pplan);
    SPI_finish();
    pfree(tr[0]);
    pfree(tr[1]);
    pfree(tr[2]);
    PG_RETURN_POINTER(result);
  }
}*/

/// RANGE QUERY
/// IN lx,ly,hx,hy,startid
PG_FUNCTION_INFO_V1(range_query);
Datum range_query(PG_FUNCTION_ARGS) 
{
  float8 bboxmin[2];
  float8 bboxmax[2];
  float8 a[2];
  float8 b[2];
  int startid,patchid;
  Vstar **tr;
  void *pplan2;
  //int res;
  int i, j, curid;
  Vstar *tmp, *tmp1;
  STACK *mystack;
  Multistar **multistar;
  
  mystack = NULL;
  tmp = palloc(sizeof(Vstar));

  SPI_connect();
  pplan2 = prepare_plan2();

  bboxmin[0] = PG_GETARG_FLOAT8(0);
  bboxmin[1] = PG_GETARG_FLOAT8(1);
  bboxmax[0] = PG_GETARG_FLOAT8(2);
  bboxmax[1] = PG_GETARG_FLOAT8(3);
  patchid = PG_GETARG_INT32(5);
  startid = PG_GETARG_INT32(4); 
  
  multistar = (Multistar**) palloc(sizeof(Multistar*)*3); // Complete fucking magic
  multistar[0] = (Multistar*) palloc(sizeof(Multistar)); // Actually storing list of 
  multistar[1] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
  multistar[2] = (Multistar*) palloc(sizeof(Multistar)); // mem locations to stars.
  multistar[0] = fetch_row(pplan2,patchid);
  elog(INFO,"RQ Got row");

  /// JUMP AND WALK HERE, TAKE CARE!
  tr = jump_and_walk(bboxmin, pplan2, startid, -1, patchid);
  SPI_finish();

  //pplan = prepare_plan(); 
  tr[0]->flag = false;
  tr[1]->flag = false;
  tr[2]->flag = false;

  elog(INFO,"RQ Jumped");


  /* bottom edge */
  a[0] = bboxmin[0];
  a[1] = bboxmin[1];
  b[0] = bboxmax[0];
  b[1] = bboxmin[1];
  walk_straight(tr, multistar, a, b, bboxmin, bboxmax, pplan2, patchid);
  /* right edge */
  a[0] = bboxmax[0];
  a[1] = bboxmin[1];
  b[0] = bboxmax[0];
  b[1] = bboxmax[1];
  walk_straight(tr, multistar, a, b, bboxmin, bboxmax, pplan2, patchid);
  /* top edge */
  a[0] = bboxmax[0];
  a[1] = bboxmax[1];
  b[0] = bboxmin[0];
  b[1] = bboxmax[1];
  walk_straight(tr, multistar, a, b, bboxmin, bboxmax, pplan2, patchid);
  /* left edge*/
  a[0] = bboxmin[0];
  a[1] = bboxmax[1];
  b[0] = bboxmin[0];
  b[1] = bboxmin[1];
  walk_straight(tr, multistar, a, b, bboxmin, bboxmax, pplan2, patchid);

  elog(INFO,"RQ Found walks");

  //multistar[0] = fetch_row(pplan2,patchid);



  //elog(INFO,"RQ Now hashing around.");
  /* now do a BFS to find all the points */
  curid = htable->sid; //latest item?
  while (curid != -1) {
    //elog(INFO,"Curid %i",curid);
    HASH_FIND_INT(htable, &curid, tmp);
    //if (tmp == NULL) {
    //elog(INFO,"Found curid! %i",tmp->sid);
    if (tmp->flag == 0) {
      //elog(INFO,"Ok...");
      tmp->flag = true;
      //elog(INFO,"Still ok!");
      for (i = 0; i < *(tmp->dim); i++) {
        //elog(INFO,"Getting %i",i);
        j = tmp->star[i];
        if (j != 0 && j < 65536) { 
          j = stitch(tmp->patchid,j); 
        }
        //elog(INFO,"Looking for id %i",j);
        if (j != 0) {
          HASH_FIND_INT(htable, &j, tmp1);
          if (tmp1 == NULL) {
            tmp1 = fetch_row_vstar2(j, multistar, pplan2, patchid);
            //tmp1 = fetch_row_vstar(j, pplan, patchid);
            //elog(INFO,"Added %i",tmp1->sid);
            HASH_ADD_INT(htable, sid, tmp1);
            //if (tmp1->sid != tmp2->sid) {
            // //elog(INFO,"Added %i, %i",tmp1->sid,tmp2->sid);
            //}
          }
          if ( (tmp1->flag == false) && (inside_rectangle(tmp1->coords, bboxmin, bboxmax) == true) ) {
            push(&mystack, j);
            //elog(INFO,"Pushed to mystack");
          }
        }
      }
    }
    //}
    //elog(INFO,"Got to popping.");
    curid = pop(&mystack);
    //elog(INFO,"Popped");
  }

  //elog(WARNING, "count %d", HASH_COUNT(htable));
  i = 0;
  for(tmp1 = htable; tmp1 != NULL; tmp1 = tmp1->hh.next) {
    if (inside_rectangle(tmp1->coords, bboxmin, bboxmax) == true) {
      //elog(WARNING, "INSIDE: %d", tmp1->id);
      i++;
      }
  }
  //elog(WARNING, "total # inside: %d", i);


  //while(htable) {
  //    tmp1 = htable;          /* copy pointer to first item     */
  //  HASH_DEL(htable, tmp1);  /* delete; users advances to next */
  //  free(tmp1);            /* optional- if you want to free  */
  //}
  /*
  HASH_ITER(hh, htable, tmp1, tmp) {
    HASH_DEL(htable,tmp1); 
   //elog(INFO,"Free %i",tmp1->sid);
    if (tmp1 != NULL) {
      free(tmp1);            
    }
  }
  */
  SPI_freeplan(pplan2);
  SPI_finish();
 //elog(INFO,"Total is %i",i);
  PG_RETURN_BOOL(true);
}


/// GET PREVIOUS IN STAR
/// QUITE SELF-EXPLAINING
int find_previous_in_star(int star[], int dim, int v)
{
  int i;
  int pos = 0;
  for (i = 0; i < dim; i++) {
    if (star[i] == v) {
      pos = i;
      break;
    }
  }
  if (pos - 1 == -1)
    return (star[dim-1]);
  else
    return (star[pos-1]);
}

int find_next_in_star(int star[], int dim, int v)
{
  int i;
  int pos = 0;
  for (i = dim; i > 0; i--) {
    if (star[i] == v) {
      pos = i;
      break;
    }
  }
  if (pos + 1 == dim)
    return (star[0]);
  else
    return (star[pos+1]);
}

/// CALCULATE TRIANGLE AREA
/// STANDARD FORMULA
double area_triangle(float8 *a, float8 *b, float8 *c)
{
  double re = ((a[0]-c[0]) * (b[1] - c[1])) - ((a[1] - c[1]) * (b[0] - c[0]));
  return (re/2);
}

/// CHECK IF LINES INTERSECT, CHECK FOUR LEFT OR RIGHT OF LINE
bool intersect(float8 *a, float8 *b, float8 *c, float8 *d)
{
  bool re = false;
  if (orient2d(a, b, c) != orient2d(a, b, d)) {
    if (orient2d(c, d, a) != orient2d(c, d, b)) {
      re = true;
    }
  }
  return re;
}

/// INTERSECTION IN 3D
/// PROBABLY TRIANGLES OR SO
void get_intersection_z(float8 *a, float8 *b, float8 *c, float8 *d, double *x)
{
	double nom, denom, u;
	nom = ((b[0] - a[0])*(a[1] - c[1])) - ((b[1] - a[1])*(a[0] - c[0]));
	denom = ((d[1] - c[1])*(b[0] - a[0])) - ((d[0] - c[0])*(b[1] - a[1]));
	u = nom / denom;
	x[0] = c[0] + (u*(d[0] - c[0]));
	x[1] = c[1] + (u*(d[1] - c[1]));
	x[2] = c[2] + (u*(d[2] - c[2]));
}

/// POINT IN BBOX
bool inside_rectangle(float8 *p, float8 *bboxmin, float8 *bboxmax)
{
bool re = false;
  if ( (p[0] >= bboxmin[0]) && (p[0] <= bboxmax[0]) ) {
    if ( (p[1] >= bboxmin[1]) && (p[1] <= bboxmax[1]) ) {
      re = true;
    }
  }
  return re;
}

/// ORIENT CHECK
/// LEFT OR RIGHT OF LINE
/// RETURNS 1 or 0
int orient2d(float8 *a, float8 *b, float8 *c)
{
  ////elog(INFO,"AB with C: %f %f %f %f and %f %f",a[0],a[1],b[0],b[1],c[0],c[1]);
  double re = ((a[0]-c[0]) * (b[1] - c[1])) - ((a[1] - c[1]) * (b[0] - c[0]));
  ////elog(INFO, "Orienting");
  if (re > 0)
    return 1; //1
  else if (re == 0) {
    ////elog(WARNING, "orient2d returns 0");
    return 0; //
  }
  else
    return -1; // -1
}


///*************************************************
///SPECIFIC TRIANGLE TYPE FUNCTIONS

PG_FUNCTION_INFO_V1(triangle_in);
Datum triangle_in(PG_FUNCTION_ARGS)
{
  char    *str = PG_GETARG_CSTRING(0);
  int    	v0, v1, v2;
  Triangle    *result;

///  //elog(WARNING, "ASCII IN");
  if (sscanf(str, " ( %d , %d , %d )", &v0, &v1, &v2) != 3)
    ereport(ERROR,
        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
         errmsg("invalid input syntax for complex: \"%s\"",
           str)));

  result = (Triangle*) palloc(sizeof(Triangle));
  result->v0 = v0;
  result->v1 = v1;
  result->v2 = v2;
  PG_RETURN_POINTER(result);
}


PG_FUNCTION_INFO_V1(triangle_out);
Datum triangle_out(PG_FUNCTION_ARGS)
{
  Triangle	*tr = (Triangle*) PG_GETARG_POINTER(0);
  char     	*result;

///  //elog(WARNING, "ASCII OUT");
  result = (char *) palloc(100);
  snprintf(result, 100, "(%d,%d,%d)", tr->v0, tr->v1, tr->v2);
  PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(triangle_recv);
Datum triangle_recv(PG_FUNCTION_ARGS)
{
  StringInfo  buf = (StringInfo) PG_GETARG_POINTER(0);
  Triangle    *result;

  //elog(WARNING, "BINARY IN");
  result = (Triangle *) palloc(sizeof(Triangle));
  result->v0 = pq_getmsgint64(buf);
  result->v1 = pq_getmsgint64(buf);
  result->v2 = pq_getmsgint64(buf);
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(triangle_send);
Datum triangle_send(PG_FUNCTION_ARGS)
{
  Triangle *tr = (Triangle *) PG_GETARG_POINTER(0);
  StringInfoData buf;

  //elog(WARNING, "BINARY OUT");
  pq_begintypsend(&buf);
  pq_sendint64(&buf, tr->v0);
  pq_sendint64(&buf, tr->v1);
  pq_sendint64(&buf, tr->v2);
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

Vec* normal(Point* p[]) {
  //float a[3],b[3];
  Vec *normal,*a,*b;
  a = palloc(sizeof(Vec));
  b = palloc(sizeof(Vec));
  normal = palloc(sizeof(Vec));

  /* Two vectors from random node of triangle */
  a->i = p[1]->x-p[0]->x;
  a->j = p[1]->y-p[0]->y;
  a->k = p[1]->z-p[0]->z;
  b->i = p[2]->x-p[0]->x;
  b->j = p[2]->y-p[0]->y;
  b->k = p[2]->z-p[0]->z;

  /* Cross product of a & b is the normal */
  normal->i = (a->j*b->k) - (b->j*a->k);
  normal->j = (a->i*b->k) - (b->i*a->k);
  normal->k = abs((a->i*b->j) - (b->i*a->j));

  return normal;
}

float slope(Vec *normal) {
  /* Returns slope in degrees */
  float s;
  float r;

  r = sqrt(pow(normal->i,2)+pow(normal->j,2));
  if (r == 0) {
    //elog(INFO,"Slope be 0");
    return 0.0;
  }
  else {
    s = atan((normal->k/r)) / PI * 180;
    //elog(INFO,"Slope be %f",s);
    return s;
  }
}

float aspect(Vec *normal) {
  if (normal->i == 0) {
    return 0.0;
  }
  else {
    return atan2(normal->j,normal->i) / PI * 180;
  }
}
/*
Create vectors from triangle or 3d plane
Vector A and B coming from B.
A = a - b 
B = c - b




float dotProduct(Vstar a, Vstar b)
{
  return a.i*b.i + a.j*b.j + a.k*b.k;
}
*/