# Integral of the magnetization along the beam: ray tracing

This note describes the former algorithm used by pathIntegral to compute, for each node of the
sensor grid, the integrals of the magnetization along the beam: a ray tracing in a CGAL AABB
tree. It is implemented in `src/pathIntegral.h` and `src/pathIntegral.cpp` (class
`PathInt::pathIntegral`), and selected with the option `--raytracing`:

```bash
pathIntegral --raytracing settings.yml
```

The default algorithm, a scan of the tetrahedrons, gives the same results 2.4 to 4 times faster;
it is described in [magnetization_integral_tetscan.md](magnetization_integral_tetscan.md). The
quantities computed (path length, integral of the magnetization, XMCD contrast) and the
preliminary steps (rotations, sensor grid) are the same for both algorithms, and are described in
sections 1 and 2 of that note.

## 1. Principle

For each node (x, y) of the sensor grid, a ray is thrown along z through the mesh. It crosses the
faces of the tetrahedrons at a sequence of points, sorted along the ray. Between two consecutive
points, the ray is inside a single tetrahedron, or outside the material. Each such segment is
integrated with the trapezoidal rule, which is exact since the magnetization is linear in each
tetrahedron (P1 finite elements).

## 2. Preparation (constructor of `pathIntegral`)

1. **Faces.** The four faces of each tetrahedron are listed as sorted triplets of node indices,
   then sorted and deduplicated: a face shared by two tetrahedrons appears once. This gives the
   vector `facets`, and the corresponding CGAL triangles (`triangles`).
2. **Owners of each face.** `facet_tets[f]` holds the (at most two) tetrahedrons having the face
   f, -1 if none: one owner for a face on the boundary of the sample, two for an inner face.
3. **AABB tree.** A CGAL AABB tree of the triangles is built once, before the parallel loop
   over the nodes, so that the tree is not built lazily by the first query inside that loop.

## 3. Intersections of a ray with the faces (`search_intersections()`)

1. **Ray.** For the node (x, y), the ray starts at (x, y, cz − 2 lz) and goes through
   (x, y, cz + 2 lz), where cz and lz are the center and the height of the bounding box of the
   mesh: it starts below the sample. Let **d** be the vector between these two points.
2. **Intersections.** `all_intersections` of the AABB tree returns all the faces crossed by the
   ray. Only the intersections that are points are kept, as hits (point P, face, abscissa
   t = **d** · (P − origin)). An intersection that is a segment (ray lying in a face) is ignored:
   the faces around it give its end points.
3. **Sort** the hits by t, then by face index.
4. **Clusters.** When the ray goes through an edge or a vertex of the mesh, all the faces sharing
   it give the same point, up to rounding errors (about 1e-22 m). Hits closer than
   `CLUSTER_TOL` × |**d**|² in t, i.e. closer than `CLUSTER_TOL` times the ray length, form a
   cluster: a single intersection point, which keeps **all** its faces.

The node is flagged inside (`flag_inside`) if the ray has at least one intersection.

## 4. Tetrahedron of a segment (`locate_tet()`)

Consider two consecutive clusters A and B, at the points PA and PB (the first hit of each
cluster). The segment [PA, PB] is inside one tetrahedron at most. The candidates are the owners
of all the faces of A and of B. For each candidate, the barycentric coordinates of the middle of
the segment are computed from signed volumes:

    λ0 = V(M, P1, P2, P3) / V,   λ1 = V(P0, M, P2, P3) / V,   ...

where V(a, b, c, d) = ((b − a) × (c − a)) · (d − a) is six times the signed volume. The candidate
with the largest smallest coordinate is chosen, if this coordinate is larger than −`BARY_TOL`;
otherwise the segment is outside the material (outside the sample, or in a cavity).

Keeping all the faces of a cluster makes this robust when the ray goes through an edge or a
vertex: the tetrahedron entered after the cluster always owns one of its faces. Before this was
done (until pull request #11), only one face was kept per point, and the tetrahedron was deduced
from the pair of consecutive faces: some segments were wrongly counted outside the material,
with errors up to 10 % on `Mx_integral` on some nodes.

## 5. Interpolation of the magnetization at a cluster point (`interpolate_u()`)

**u** is interpolated on the first face of the cluster, with the area coordinates of the point
P in the triangle (P0, P1, P2):

    αi = Ai / (A0 + A1 + A2),   **u**(P) = Σ αi **u**i

where Ai is the area of the triangle formed by P and the two vertices other than Pi. On a face,
this is the restriction of the P1 interpolation of each tetrahedron sharing it, so the value does
not depend on the face chosen within the cluster.

## 6. Integration along the ray (`path_integrals()`)

For each segment between consecutive clusters, of length δ = |PB − PA|, with **u**A and **u**B
interpolated at its ends, and inside = 1 if a tetrahedron was found, 0 otherwise:

    I        += ½ (uA + uB) δ inside J / µ0
    arg+     += ½ (λ+ + λ− + (λ+ − λ−) s) δ inside
    arg−     += ½ (λ+ + λ− − (λ+ − λ−) s) δ inside
    s         = ½ k · (uA + uB)

where k is the unit vector of the ray, and J, λ+ = `lbp` and λ− = `lbm` (given in nm⁻¹, multiplied
by 10⁹) are the parameters of the volume region of the tetrahedron.

- **Path length**: if `filled` is false, L = Σ δ inside. If `filled` is true, every segment counts,
  so L is the distance from the first to the last intersection point.
- **XMCD contrast**: C = (I+ − I−) / (I+ + I−), with I± = exp(−arg±).

The helper `integrate_u()` computes the same path length and the integral of **u**, without the
material parameters; it is not used by pathIntegral.

## 7. Tolerances

| Constant | Value | Role |
|---|---|---|
| `CLUSTER_TOL` | 1e-9 | relative to the ray length: hits closer than this are the same point |
| `BARY_TOL` | 1e-6 | tolerance on the barycentric coordinates of the middle of a segment |

## 8. Cost and parallelism

The loop over the grid nodes is parallel (OpenMP, dynamic schedule by chunks of 256 nodes, in
`main.cpp`); the intersections and hits are stored in `thread_local` vectors, reused from one ray
to the next.

- A ray that misses the sample is rejected near the root of the AABB tree, in about 0.1 µs.
- A ray that crosses the sample costs about 60 to 75 µs on the test cases, for some forty faces
  crossed: the computation of the intersections and their sorting dominate.

On the sphere of `simulations/oneSphere_D100nm_Js1T_v1` (91 275 tetrahedra, 160 801 grid nodes),
the ray tracing takes about 0.7 s on 8 threads, the AABB tree included, against about 0.2 s for
the tetrahedron scan.

## 9. Differences with the tetrahedron scan

- The ray tracing does not assume that the beam is along z, nor that the grid is regular: each
  ray only depends on its node. It assumes, like the scan, a conforming mesh.
- `flag_inside` is set as soon as the ray has an intersection, so a ray tangent to the sample, with
  a zero path length, may or may not be flagged depending on rounding errors. The tetrahedron scan
  accepts such contacts up to a tolerance, so that they are always flagged.
- Otherwise, the results agree to about 1e-9 (relative).
