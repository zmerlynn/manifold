# Worsened-case fuzz seeds for RemoveSelfIntersections

317 cases discovered by `scripts/fuzz_overnight.sh 1 5000` where
the pipeline made the input WORSE (= `post_pierces > pre_pierces`).
These are deterministic reproducers — given the spike at commit
`20a177c4` or compatible, running the listed (master_seed, class,
case) tuple regenerates the exact input.

Run a single case:
```
OVERLAP3D_FUZZ_SEED=<master_seed> build/extras/overlap3d_proto --advfuzz
```
The summary at the end will show the per-class outcomes for that
master_seed; the specific case is the `case_index` within the class
loop.

## Class index

| class | name | description |
|---|---|---|
| 0 | shallow-tight | Add of cube + rotated cube, displaced to kPow=30 |
| 1 | near-coplanar-slabs | Add of two thin slabs at small tilt |
| 2 | 3-cube-chain | chained Add of three rotated cubes |
| 3 | 4-cube-chain | chained Add of four rotated cubes |
| 4 | 5-cube-chain | chained Add of five rotated cubes |
| 5 | 6-cube-chain | chained Add of six rotated cubes |
| 6 | vertex-on-face | cube + cube positioned so vert lands on face |
| 7 | subtract-rotated | Subtract of rotated cube from cube |
| 8 | subtract-3-chain | three-stage Add+Subtract chain |
| 9 | smith-three-mutual | three orthogonal slabs (Smith UCAM-CL-TR-766 fig 9.1) |
| 10 | smith-four-coplanar | four near-coplanar slabs (Smith p.84 field bug) |
| 11 | cospherical-shell | shell carved out of sphere |
| 12 | cgal-12-cube-chain | 12-deep nested Boolean of cubes (CGAL stress case) |

## Failure breakdown by class

```
class7  (subtract-rotated):   182 worsened cases  (57%)
class8  (subtract-3-chain):    52              (16%)
class12 (cgal-12-cube-chain):  37              (12%)
class9  (smith-three-mutual):  17               (5%)
class5  (6-cube-chain):        15               (5%)
class3  (4-cube-chain):         8               (3%)
class4  (5-cube-chain):         6               (2%)
class0/1/2/6/10/11:             0
```

**Subtract operations account for 74% of worsened cases.** This
matches the cray pathology (= Subtract back-side normal flip
confuses the orientation classifier) and is the highest-leverage
class to investigate.

## Worst (biggest delta) representatives

These have stable `.obj` copies in `test/models/fuzz/`:

| file | input pierces | output pierces | delta | class |
|---|---|---|---|---|
| `worst_delta_sub3chain.obj` | 2 | 9 | **+7** | subtract-3-chain (seed2227, case175) |
| `subrotated_pre2_post8.obj` | 2 | 8 | +6 | subtract-rotated (seed4816, case139) |
| `smith_three_mutual_pre4_post8.obj` | 4 | 8 | +4 | smith-three-mutual (seed1986, case111) |
| `cgal12cube_pre2_post6.obj` | 2 | 6 | +4 | cgal-12-cube-chain (seed879, case50) |
| `subrotated_pre3_post5.obj` | 3 | 5 | +2 | subtract-rotated (seed5, case13) |
| `cube4chain_pre2_post4.obj` | 2 | 4 | +2 | 4-cube-chain (seed135, case47) |

`subrotated_pre3_post5.obj` is the simplest representative — single
Subtract operation, only 3 pierces in, smallest mesh. Good starting
point for debugging.

## Full worsened case list

Format: `master_seed,class,case_index,pre_pierces,post_pierces,delta`

```csv
5,class7,13,3,5,2
19,class7,47,4,5,1
42,class7,14,4,5,1
86,class7,80,4,5,1
87,class12,71,3,6,3
91,class7,26,2,3,1
102,class7,58,2,3,1
103,class7,123,3,4,1
135,class3,47,2,4,2
135,class8,81,8,8,0
136,class12,23,3,3,0
156,class9,70,6,7,1
158,class7,15,2,3,1
183,class5,92,5,6,1
201,class9,108,2,4,2
206,class7,82,2,3,1
262,class7,81,4,5,1
286,class7,129,2,3,1
291,class7,149,3,4,1
301,class7,123,2,3,1
314,class3,129,3,4,1
317,class8,92,5,6,1
341,class8,114,2,4,2
346,class12,50,2,3,1
357,class5,148,6,7,1
361,class9,23,4,5,1
369,class7,65,3,7,4
411,class9,11,3,4,1
447,class9,165,4,5,1
447,class7,41,3,4,1
462,class7,196,2,3,1
465,class8,142,2,4,2
479,class7,128,3,5,2
540,class7,164,3,5,2
554,class7,113,2,4,2
555,class12,156,2,3,1
561,class7,33,2,5,3
574,class7,103,2,5,3
577,class3,4,4,6,2
584,class7,123,3,5,2
593,class8,57,9,11,2
604,class8,194,5,8,3
613,class12,99,3,4,1
624,class3,16,3,4,1
650,class12,50,3,4,1
653,class7,29,2,3,1
667,class7,86,3,4,1
682,class9,33,4,5,1
697,class7,68,5,6,1
700,class7,29,3,4,1
720,class5,177,2,3,1
752,class12,123,3,4,1
786,class7,118,2,3,1
786,class7,124,2,3,1
792,class7,55,4,5,1
816,class8,131,4,7,3
820,class3,148,2,3,1
822,class7,69,3,4,1
824,class7,128,3,4,1
840,class7,170,2,3,1
848,class9,46,2,3,1
865,class12,151,2,4,2
877,class12,108,2,3,1
878,class7,80,4,5,1
879,class12,50,2,6,4
887,class12,18,3,4,1
892,class12,18,3,4,1
905,class7,108,4,5,1
905,class7,159,2,3,1
909,class7,179,2,4,2
913,class7,38,2,3,1
934,class7,68,4,6,2
934,class7,69,4,5,1
1031,class12,30,3,4,1
1031,class7,189,2,4,2
1048,class3,9,3,4,1
1056,class12,156,2,4,2
1058,class7,52,3,8,5
1075,class7,28,3,4,1
1086,class7,118,4,5,1
1090,class7,28,3,4,1
1090,class8,11,3,5,2
1094,class12,86,4,5,1
1100,class7,141,3,4,1
1102,class7,68,3,4,1
1110,class7,170,3,4,1
1114,class9,156,4,5,1
1129,class8,73,2,3,1
1140,class8,131,3,4,1
1148,class8,170,2,6,4
1156,class12,98,3,4,1
1170,class7,79,3,4,1
1190,class8,11,2,5,3
1198,class8,159,5,8,3
1207,class7,99,3,4,1
1213,class7,116,2,3,1
1227,class7,21,4,5,1
1235,class12,98,2,4,2
1294,class7,103,2,3,1
1306,class7,134,3,5,2
1311,class9,82,4,6,2
1339,class7,110,3,8,5
1356,class7,8,4,5,1
1360,class12,12,2,3,1
1392,class7,98,3,4,1
1399,class7,118,2,3,1
1402,class3,29,4,5,1
1408,class12,156,2,4,2
1421,class7,69,3,4,1
1426,class7,134,3,4,1
1435,class7,67,2,5,3
1462,class7,144,2,3,1
1471,class3,75,4,5,1
1479,class7,114,3,4,1
1485,class7,72,3,4,1
1498,class7,11,5,6,1
1505,class8,86,4,7,3
1508,class9,28,2,4,2
1521,class7,82,2,3,1
1593,class7,123,3,4,1
1611,class12,156,3,5,2
1614,class7,127,2,3,1
1625,class7,176,3,4,1
1648,class7,53,3,4,1
1671,class12,103,4,5,1
1680,class12,124,4,5,1
1689,class7,29,3,4,1
1696,class9,182,2,3,1
1707,class7,170,3,4,1
1716,class7,128,3,4,1
1721,class7,176,2,3,1
1722,class7,148,3,4,1
1738,class7,113,2,3,1
1740,class9,79,5,7,2
1745,class7,67,5,6,1
1747,class7,67,2,3,1
1763,class7,68,3,4,1
1766,class12,151,2,3,1
1778,class7,134,3,4,1
1785,class7,90,2,3,1
1791,class7,80,3,4,1
1838,class7,128,3,4,1
1862,class7,29,2,3,1
1868,class9,49,3,4,1
1875,class7,11,2,3,1
1889,class12,98,2,4,2
1894,class7,87,4,5,1
1953,class7,73,2,3,1
1981,class12,84,2,3,1
1986,class9,111,4,8,4
2001,class8,96,16,15,-1
2002,class7,16,2,3,1
2019,class7,29,2,3,1
2024,class7,134,2,4,2
2032,class12,108,4,5,1
2079,class7,84,3,4,1
2092,class7,75,3,4,1
2106,class7,164,2,3,1
2115,class7,166,2,3,1
2189,class5,44,25,14,-11
2227,class8,175,2,9,7
2257,class7,131,2,3,1
2278,class12,156,3,4,1
2280,class7,127,2,3,1
2317,class9,32,2,3,1
2341,class7,80,4,5,1
2342,class7,179,2,3,1
2350,class8,21,2,6,4
2351,class7,127,3,7,4
2398,class7,33,3,4,1
2398,class7,75,2,3,1
2418,class7,11,4,5,1
2428,class7,15,2,4,2
2441,class7,176,3,4,1
2487,class7,33,2,4,2
2540,class7,145,2,8,6
2550,class7,86,3,4,1
2571,class7,40,2,3,1
2585,class7,189,4,5,1
2588,class12,156,4,5,1
2603,class7,135,5,6,1
2622,class7,3,4,5,1
2645,class7,33,2,4,2
2673,class7,83,2,6,4
2717,class7,78,3,8,5
2728,class12,99,3,4,1
2737,class7,170,3,4,1
2741,class9,17,2,4,2
2780,class12,5,4,5,1
2787,class7,128,3,4,1
2795,class7,87,3,4,1
2810,class7,82,3,4,1
2811,class3,99,2,3,1
2827,class7,87,3,4,1
2853,class7,38,2,3,1
2876,class7,93,2,3,1
2877,class3,142,4,5,1
2877,class12,124,2,4,2
2900,class7,87,2,3,1
2933,class7,108,4,5,1
2944,class7,108,3,4,1
2960,class8,127,5,8,3
2976,class7,69,2,5,3
2983,class8,75,3,4,1
2987,class7,80,3,5,2
3000,class8,176,5,8,3
3019,class7,159,2,3,1
3060,class7,148,2,3,1
3072,class8,6,5,9,4
3080,class7,160,2,3,1
3115,class7,182,2,3,1
3115,class7,53,2,3,1
3132,class7,131,3,4,1
3155,class7,87,4,5,1
3155,class7,90,3,4,1
3194,class7,18,3,4,1
3243,class7,103,3,8,5
3261,class9,114,2,3,1
3267,class9,191,4,6,2
3278,class12,71,2,3,1
3286,class12,123,2,3,1
3307,class7,38,3,7,4
3326,class7,11,2,3,1
3343,class7,114,2,3,1
3360,class7,98,3,4,1
3401,class7,21,3,4,1
3401,class7,124,2,3,1
3421,class9,83,2,4,2
3438,class7,128,3,4,1
3552,class7,33,4,6,2
3554,class7,55,3,4,1
3562,class12,0,5,8,3
3596,class7,103,3,4,1
3617,class7,165,2,3,1
3617,class7,179,3,4,1
3621,class7,103,4,5,1
3673,class7,134,2,3,1
3680,class7,84,3,7,4
3713,class5,30,7,8,1
3736,class7,128,3,4,1
3755,class12,151,3,5,2
3759,class7,80,2,3,1
3759,class9,32,3,4,1
3782,class3,75,5,6,1
3786,class12,156,3,5,2
3789,class12,156,2,4,2
3790,class7,38,2,4,2
3793,class9,72,3,4,1
3823,class7,68,2,3,1
3851,class7,49,3,4,1
3858,class8,120,6,9,3
3936,class7,53,5,6,1
3942,class7,28,2,3,1
3944,class8,160,7,8,1
3954,class7,86,4,5,1
3989,class7,118,4,5,1
4012,class7,69,4,5,1
4022,class7,134,2,3,1
4030,class7,116,2,3,1
4040,class7,28,2,3,1
4055,class7,80,3,4,1
4117,class7,148,2,3,1
4148,class7,42,5,6,1
4187,class7,16,3,4,1
4218,class7,38,3,4,1
4242,class7,176,3,4,1
4254,class7,164,2,3,1
4270,class7,170,2,3,1
4297,class7,103,3,4,1
4308,class7,164,4,5,1
4338,class7,87,2,3,1
4361,class7,170,2,3,1
4377,class7,33,4,5,1
4413,class9,32,5,6,1
4438,class7,15,2,3,1
4488,class7,33,5,6,1
4490,class7,22,3,8,5
4492,class7,131,3,4,1
4500,class7,8,3,4,1
4502,class8,49,6,9,3
4503,class3,32,3,4,1
4516,class12,151,4,5,1
4555,class7,118,2,3,1
4587,class7,98,5,6,1
4602,class7,69,3,4,1
4610,class7,67,2,3,1
4674,class8,41,5,8,3
4711,class12,86,3,5,2
4716,class8,163,18,17,-1
4729,class12,156,3,4,1
4774,class7,149,3,7,4
4774,class7,90,2,3,1
4776,class12,123,3,4,1
4816,class7,139,2,8,6
4830,class7,3,3,4,1
4834,class7,123,3,4,1
4839,class7,16,2,3,1
4839,class8,194,4,7,3
4874,class7,85,4,7,3
4885,class8,78,5,9,4
4902,class7,38,4,5,1
4915,class7,182,2,3,1
4937,class7,28,4,5,1
4970,class8,9,4,12,8
4973,class7,177,3,4,1
4984,class4,49,2,3,1
4990,class7,91,4,5,1
```

## What's already in `test/models/fuzz/`

```
worst_delta_sub3chain.obj          (seed2227, class8, case175, pre=2, post=9, +7)
subrotated_pre2_post8.obj          (seed4816, class7, case139, pre=2, post=8, +6)
smith_three_mutual_pre4_post8.obj  (seed1986, class9, case111, pre=4, post=8, +4)
cgal12cube_pre2_post6.obj          (seed879,  class12, case50, pre=2, post=6, +4)
subrotated_pre3_post5.obj          (seed5,    class7, case13,  pre=3, post=5, +2)
cube4chain_pre2_post4.obj          (seed135,  class3, case47,  pre=2, post=4, +2)
```

## Debugging entry points

The pipeline's **pierce-monotonicity gate** in
`src/overlap_removal.cpp::RunOverlapRemovalImpl` is supposed to
detect post > pre and fall back to input. The 317 worsened cases
are precisely cases where this gate failed to catch the regression.

Two possible causes:
1. The gate compares pierce counts on the FINAL output. If the
   pipeline produces a NoError manifold whose pierce count exceeds
   the input, the gate should fall back. Verify by running the gate
   logic on these cases.
2. The pierce-counting predicate (`CheckSelfIntersection`) uses a
   relTol-based test; some "real" pierces might be missed in the
   input but caught in the output (= apparent regression). Less
   likely but possible.

`subrotated_pre3_post5.obj` is the most tractable starting case:
single Subtract operation, ~hundred tris, only 3 pierces in. Load
it, run RemoveSelfIntersections, observe what comes out, trace
why the gate let it through.
