

CREATE OR REPLACE FUNCTION npoints(integer[])
  RETURNS integer AS
'example.so', 'npoints'
  LANGUAGE c IMMUTABLE STRICT;
  
  CREATE OR REPLACE FUNCTION getstar(integer,integer[], bytea)
  RETURNS integer[] AS
'example.so', 'getstar'
  LANGUAGE c IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION getx(integer,integer[],bytea)
  RETURNS double precision AS
'example.so', 'getpoint'
  LANGUAGE c IMMUTABLE STRICT;


CREATE OR REPLACE FUNCTION getpoint(IN integer,IN integer[],IN bytea,
    OUT x real, OUT y real, OUT z real)
    AS 'MODULE_PATHNAME', 'getpoint'
    LANGUAGE C IMMUTABLE STRICT;
    
CREATE OR REPLACE FUNCTION getrow(IN integer,IN integer[],IN bytea,IN bytea,
    OUT id integer, OUT x double precision, OUT y double precision, OUT z double precision, OUT stars integer[])
    AS 'MODULE_PATHNAME', 'getrow'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION getpoints(IN integer[], IN bytea,
	OUT id integer, OUT x real, OUT y real, OUT z real)
  RETURNS SETOF record 
  AS 'MODULE_PATHNAME', 'getpoints'
  LANGUAGE c IMMUTABLE STRICT  ;

CREATE OR REPLACE FUNCTION simplify(IN integer[], IN bytea, IN bytea,
	OUT x double precision, OUT y double precision, OUT z double precision, OUT dif double precision)
  RETURNS SETOF record 
  AS 'MODULE_PATHNAME', 'simplify'
  LANGUAGE c IMMUTABLE STRICT  ;
  

CREATE OR REPLACE FUNCTION gettin(IN integer[], IN bytea, IN bytea,
	OUT id integer, OUT x real, OUT y real, OUT z real, OUT star integer[])
  RETURNS SETOF record 
  AS 'MODULE_PATHNAME', 'gettin'
  LANGUAGE c IMMUTABLE STRICT  ;
  
CREATE OR REPLACE FUNCTION gettin2(IN integer[], IN bytea, IN bytea,
	OUT id integer, OUT x double precision, OUT y double precision, OUT z double precision, OUT star integer[])
  RETURNS SETOF record 
  AS 'MODULE_PATHNAME', 'gettin'
  LANGUAGE c IMMUTABLE STRICT  ;
  
  CREATE OR REPLACE FUNCTION gettin_convex(IN integer[], IN bytea, IN bytea,
	OUT id integer, OUT x double precision, OUT y double precision, OUT z double precision, OUT star integer[])
  RETURNS SETOF record 
  AS 'MODULE_PATHNAME', 'gettin_convex'
  LANGUAGE c IMMUTABLE STRICT  ;
  
CREATE TYPE triangle;

CREATE OR REPLACE FUNCTION triangle_in (cstring)
  RETURNS triangle
AS 'MODULE_PATHNAME', 'triangle_in'
  LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION triangle_out (triangle)
  RETURNS cstring
AS 'MODULE_PATHNAME', 'triangle_out'
  LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION triangle_recv(internal)
  RETURNS triangle
AS 'MODULE_PATHNAME', 'triangle_recv'
  LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION triangle_send(triangle)
  RETURNS bytea
AS 'MODULE_PATHNAME', 'triangle_send'
  LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE triangle 
(
 internallength = 12,
 input = triangle_in,
 output = triangle_out,
 receive = triangle_recv,
 send = triangle_send
);

CREATE OR REPLACE FUNCTION testmulti(integer,integer) RETURNS integer
  AS 'MODULE_PATHNAME', 'testmulti'
  LANGUAGE c IMMUTABLE STRICT;

  
CREATE OR REPLACE FUNCTION pl2(double precision, double precision, integer)
  RETURNS bytea AS
'MODULE_PATHNAME', 'pl'
  LANGUAGE c IMMUTABLE STRICT
  
CREATE OR REPLACE FUNCTION profile_count_intersections(double precision, double precision, double precision, double precision, integer)
    RETURNS integer
    AS 'MODULE_PATHNAME', 'profile_count_intersections'
    LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION range_query(double precision, double precision, double precision, double precision, integer, integer)
    RETURNS boolean 
    AS 'MODULE_PATHNAME', 'range_query'
    LANGUAGE C STRICT IMMUTABLE;
