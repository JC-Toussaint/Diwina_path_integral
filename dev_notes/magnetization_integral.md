# Integral of the magnetization along the beam

This note describes how pathIntegral computes, for each node of the sensor grid, the integrals
of the magnetization along the beam. The algorithm is implemented in `src/tetScan.h` and
`src/tetScan.cpp` (class `PathInt::tetScan`), and called from `src/main.cpp`.

## 1. Problem

The micromagnetic sample is a tetrahedral mesh. The magnetization is known at the mesh nodes, as
the reduced magnetization **u** (|**u**| = 1), and is linear in each tetrahedron (P1 finite
elements). The magnetization of a tetrahedron is **M** = (J/µ0) **u**, where J = µ0 Ms is the
`Js` of its volume region.

After the rotations, the beam is along z. For each node (x, y) of the sensor grid, the vertical
line through (x, y) crosses the sample, and pathIntegral computes:

- the path length L(x, y), the length of the line inside the material;
- the integral of the magnetization, **I**(x, y) = ∫ **M**(x, y, z) dz, in A;
- the XMCD contrast C(x, y), from the absorption of the two circular polarizations.

These are the columns `path_length`, `Mx_integral`, `My_integral`, `Mz_integral` and `contrast`
of `*_M_integrals.out`.

## 2. Preliminary steps (`main.cpp`)

1. **Mesh and magnetization.** The node coordinates are multiplied by `length_unit`, so all the
   lengths are in meters. **u** is read from the `.sol` file, or computed from the expressions
   of `initial_magnetization`.
2. **Rotations.** The two rotations of the settings (`angle1` around `axe1`, then `angle2` around
   `axe2`) are applied to the node positions **and** to **u**. The bounding box of the mesh
   (center c, lengths l) is then updated.
3. **Sensor grid** (`Fem2d`). A regular square grid of N × N nodes, centered on the projection
   of the sample on the xy plane. Its side is max(lx, ly) × `relative_size`, and its step is
   `pixel_size`. Node (ix, iy) is at (xmin + ix·p, ymin + iy·p), and has index ix·N + iy.

## 3. Principle: scan the tetrahedrons

All the integration lines are vertical, so instead of throwing one ray per node through the mesh,
the algorithm works tetrahedron by tetrahedron:

- each tetrahedron is projected on the xy plane: only the grid nodes inside its projection can
  see it;
- for such a node, the vertical line crosses the tetrahedron along a segment, the **chord**
  [z_in, z_out], computed from the barycentric coordinates;
- since **u** is linear in the tetrahedron, it is linear along the chord, and the trapezoidal
  rule on the chord is exact.

The integral along a line is then the sum over the chords of the line.

## 4. Barycentric coordinates

For a tetrahedron with vertices P0, P1, P2, P3, let J = [P1 − P0, P2 − P0, P3 − P0] be the 3 × 3
matrix of its edges from P0. The barycentric coordinates of a point X are

    (λ1, λ2, λ3) = J⁻¹ (X − P0),    λ0 = 1 − λ1 − λ2 − λ3

so each λi is an affine function of X:

    λi(X) = δi0 + gi · (X − P0)

where g1, g2, g3 are the rows of J⁻¹, and g0 = −(g1 + g2 + g3). X is inside the tetrahedron if and
only if λi(X) ≥ 0 for i = 0..3, and **u**(X) = Σ λi(X) **u**i.

The constructor of `tetScan` computes, for each tetrahedron and in parallel: P0, the four
gradients gi (`TetData::grad`, a 4 × 3 matrix), and the bounding box of its projection on the xy
plane (`xmin`, `xmax`, `ymin`, `ymax`).

## 5. Candidate tetrahedrons of each node

A node can only be crossed by the tetrahedrons whose projected bounding box contains it. For
each tetrahedron, the grid indices covered by its bounding box are

    ix0 = max(0, floor((xmin − x0)/p)),   ix1 = min(N − 1, ceil((xmax − x0)/p))

and similarly for iy, where (x0, y0) is the first grid node and p the pixel size. The rounding
is outwards, so no node is missed; the extra candidates are rejected by the chord computation.

The candidates are stored per node in compressed sparse row format: `first[n]` is the position
of the first candidate of node n in the array `candidates`, and `first[n+1]` the end. The
construction is parallel, without atomic operations, by giving each grid row ix to a single
thread:

1. **Tetrahedrons of each row**: for each tetrahedron, compute its index ranges, and count it in
   each row ix0..ix1, with one counter array per thread. The counters are turned into positions,
   and a second pass over the tetrahedrons, with the same static distribution over the threads,
   fills the list of the tetrahedrons of each row.
2. **Number of candidates of each node**: for each row, in parallel, count the candidates of its
   nodes from the tetrahedrons of the row; a prefix sum gives `first`.
3. **Candidates**: for each row, in parallel, write the candidates of its nodes.

For each node, the candidates are sorted by increasing tetrahedron index, whatever the number of
threads.

## 6. Chord of a vertical line in a tetrahedron (`chord()`)

Along the vertical line through (x, y), write z = P0z + w. Each barycentric coordinate is an
affine function of w:

    λi = ai + bi w,   ai = δi0 + gix (x − P0x) + giy (y − P0y),   bi = giz

The line is inside the tetrahedron where all the λi are non-negative:

- if bi > 0, λi ≥ 0 gives w ≥ −ai/bi: lower bound;
- if bi < 0, λi ≥ 0 gives w ≤ −ai/bi: upper bound;
- if bi ≈ 0 (|bi| ≤ `GRAD_TOL` |gi|), the line is parallel to the face opposite to vertex i:
  λi is constant along the line, and the line misses the tetrahedron if ai < −`BARY_TOL`.

So w_in = max of the lower bounds, w_out = min of the upper bounds, and the chord is
[P0z + w_in, P0z + w_out]. The line misses the tetrahedron if w_in > w_out + tol_z.

A line through an edge or a vertex of the tetrahedron gives w_in = w_out up to rounding errors.
The tolerance tol_z makes such a contact count as a zero-length chord, so the result does not
depend on rounding: w_out is then set to max(w_in, w_out).

## 7. Integration along a line

For each node n, in parallel over the nodes:

1. **Chords**: compute the chord of the line in each candidate tetrahedron. The node is flagged
   inside (`flag_inside`) if at least one candidate is touched, even with a zero-length chord.
   Only the chords longer than tol_z are kept.
2. **Sort** the chords by (z_in, z_out, tetrahedron index).
3. **Duplicates**: a chord starting before the end of the previous kept chord (z_in < z_last −
   tol_z) is dropped. This happens when the line goes along a vertical face shared by two
   tetrahedrons: both give the same chord. In a conforming mesh the kept chords then cover the
   line inside the material exactly once.
4. **Integrals**: for each kept chord of tetrahedron t, with δ = z_out − z_in,
   **u**in = **u**(x, y, z_in) and **u**out = **u**(x, y, z_out) interpolated in t:

       I        += ½ (uin + uout) δ J / µ0
       arg+     += ½ (λ+ + λ− + (λ+ − λ−) s) δ
       arg−     += ½ (λ+ + λ− − (λ+ − λ−) s) δ
       s         = ½ (uin,z + uout,z)

   where J, λ+ = `lbp` and λ− = `lbm` are the parameters of the volume region of t. `lbp` and
   `lbm` are given in nm⁻¹ in the settings and multiplied by 10⁹. s is the mean of u·k along the
   chord, k = z being the beam direction.
5. **Path length**: if `filled` is false, L is the sum of the δ (material only). If `filled` is
   true, the hollow parts crossed by the line count as material: L = z_last − z_first, from the
   entry in the first chord to the exit of the last one.
6. **XMCD contrast**: with the transmitted intensities I± = exp(−arg±) of the two circular
   polarizations,

       C = (I+ − I−) / (I+ + I−)

All the integrands are linear along a chord, so the trapezoidal rule is exact for the P1
magnetization.

## 8. Tolerances

| Constant | Value | Role |
|---|---|---|
| `Z_TOL` | 1e-9 | tol_z = `Z_TOL` × lz (height of the mesh): zero-length chords, contacts, overlapping chords |
| `GRAD_TOL` | 1e-12 | relative threshold on giz under which the line is parallel to a face |
| `BARY_TOL` | 1e-10 | tolerance on λi for a line parallel to a face |

## 9. Cost and parallelism

- The preparation (barycentric gradients, candidate lists) costs O(number of tetrahedrons ×
  pixels covered by a tetrahedron), and is parallel without atomic operations.
- The integration costs O(number of candidates), with a sort of the chords of each line. It is
  parallel over the nodes, with a dynamic schedule since the nodes above the sample cost more
  than the others.
- Nodes outside the projection of the sample have no candidate, and cost almost nothing.

On the sphere of `simulations/oneSphere_D100nm_Js1T_v1` (91 275 tetrahedra, 160 801 grid nodes),
the whole step takes about 0.2 s on 8 threads, against about 0.7 s for the former ray tracing.

## 10. Assumptions

- The beam is along z after the rotations.
- The sensor grid is regular and square, with node (ix, iy) at index ix·N + iy.
- The mesh is conforming: two neighbouring tetrahedrons share a whole face. The rule that drops
  the duplicate chords relies on it.
- The magnetization is P1, and the material parameters (Js, lbp, lbm) are constant in each volume
  region.

## 11. Validation

- With a uniform magnetization, `Mx_integral = path_length × Js / µ0` holds on every node to
  about 3e-11 (relative), including the lines through edges and vertices of the mesh.
- On a hollow sphere (radii 100 and 200 nm, 10 nm mesh), the path length agrees with the analytic
  chord to within 0.6 nm, in the cavity as well as in the shell, the remaining error coming from
  the faceting of the spheres.
- The results agree with the former ray tracing algorithm to about 1e-9 (relative).

## 12. Former algorithm: ray tracing

The former algorithm is still available with `pathIntegral --raytracing settings.yml`
(`src/pathIntegral.h`, `src/pathIntegral.cpp`). For each grid node, it throws a ray through a CGAL
AABB tree of all the faces of the tetrahedrons, sorts the intersection points along the ray,
groups those closer than 1e-9 times the ray length (ray through an edge or a vertex), and finds
the tetrahedron of each segment as the one, among the tetrahedrons sharing the faces at its ends,
that contains its middle. It gives the same results, but is 2.4 to 4 times slower.

## 13. What comes next

The integrals of the magnetization are the input of the holographic phase. In `main.cpp`, the
in-plane integrated magnetization is rotated by 90° around z (M′ = ez × M, i.e. Mx′ = −My,
My′ = Mx), and the magnetic part of the phase is computed by a 2D fast multipole method
(`pot2D::scalfmm2d_sum`, `src/scalfmm2d_sum.cpp`) over the sensor grid.
