# multistar 0.1 Alpha
*Here be dragons*
Multistar data structure for storing TINs in PostgreSQL packaged as an extension. Currently has some issues, but works on Point and Line Intersection.

## Data used
See the thesis_tools repository for the file multistar.py, which can create the data structure queried. At the moment the tablename is hardcoded once in the source code, which should be changed according to the name used in the generator python file. Although the table construction is done by Python, the command used is:

```SQL
CREATE TABLE multistar (id int, bbox box3d, offsets int[], points bytea, stars bytea);
ALTER TABLE multistar ADD PRIMARY KEY (id);
CREATE INDEX idx_multistar ON {} using GIST(ST_FORCE_3DZ(bbox));
CREATE OR REPLACE VIEW multistar_all AS SELECT id,st_force_3dz(bbox::geometry) AS st_force_3dz
FROM multistar;
```


## Use
After compiling with the provided Makefile and loading the provided extension, it is used as follows for a profile, returning a LineStringZ geometry.

```SQL
SELECT profile_count_intersections
(82348.405,449248.327,82556.899,449248.503,id) FROM multistar_l6
WHERE bbox && ST_GEOMFROMTEXT( ’POINT(82348.405 ,449248.327) ’ )
```

It can also be used for point location, returning the TriangleZ in which the requested points is located.

```SQL
SELECT pl(82419.172,449005.397,id) FROM multistar_l4 WHERE bbox &&
ST_POINT(82419.172 449005.397);
```