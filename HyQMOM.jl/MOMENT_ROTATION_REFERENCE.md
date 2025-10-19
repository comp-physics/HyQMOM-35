# Moment Rotation Reference for Z-axis (+90°)

## Coordinate Transformation
```
(x, y, z) → (y, -x, z)
(u, v, w) → (v, -u, w)
```

## Moment Transformation Rule
```
M_{ijk}^new = ⟨u'^i v'^j w'^k⟩
            = ⟨v^i (-u)^j w^k⟩
            = (-1)^j ⟨v^i u^j w^k⟩
            = (-1)^j M_{jik}^old
```

## Full 35-Moment Mapping

```
Index | Old Moment | New Moment | Transform              | Notes
------|------------|------------|------------------------|-------
  1   | M000       | M000       | M000                   | scalar, invariant
  2   | M100       | M010       | -M010 (= (-1)^0 M010)  | u → v
  3   | M200       | M020       | M020                   | u² → v²
  4   | M300       | M030       | -M030                  | u³ → -v³
  5   | M400       | M040       | M040                   | u⁴ → v⁴
  6   | M010       | M100       | -M100                  | v → -u
  7   | M110       | M110       | -M110                  | uv → -v(-u) = -uv
  8   | M210       | M120       | M120                   | u²v → v²(-u) = -uv²
  9   | M310       | M130       | -M130                  | u³v → v³(-u) = -u³v
 10   | M020       | M200       | M200                   | v² → u²
 11   | M120       | M210       | -M210                  | uv² → v(-u)² = vu²
 12   | M220       | M220       | M220                   | u²v² → v²u²
 13   | M030       | M300       | -M300                  | v³ → -u³
 14   | M130       | M310       | M310                   | uv³ → v(-u)³ = uv³
 15   | M040       | M400       | M400                   | v⁴ → u⁴
 16   | M001       | M001       | M001                   | w invariant
 17   | M101       | M011       | -M011                  | uw → vw
 18   | M201       | M021       | M021                   | u²w → v²w
 19   | M301       | M031       | -M031                  | u³w → -v³w
 20   | M002       | M002       | M002                   | w² invariant
 21   | M102       | M012       | -M012                  | uw² → vw²
 22   | M202       | M022       | M022                   | u²w² → v²w²
 23   | M003       | M003       | M003                   | w³ invariant
 24   | M103       | M013       | -M013                  | uw³ → vw³
 25   | M004       | M004       | M004                   | w⁴ invariant
 26   | M011       | M101       | -M101                  | vw → -uw
 27   | M111       | M111       | M111                   | uvw → v(-u)w = -uvw, but signs...
 28   | M211       | M121       | -M121                  | u²vw → v²(-u)w = -v²uw
 29   | M021       | M201       | M201                   | v²w → u²w
 30   | M121       | M211       | -M211                  | uv²w → v(-u)²w = vu²w
 31   | M031       | M301       | -M301                  | v³w → -u³w
 32   | M012       | M102       | -M102                  | vw² → -uw²
 33   | M112       | M112       | M112                   | uvw² → v(-u)w² = -uvw²
 34   | M013       | M103       | -M103                  | vw³ → -uw³
 35   | M022       | M202       | M202                   | v²w² → u²w²
```

## Key Patterns
1. **Invariant under rotation (no sign change):**
   - Pure w moments (w^k): M001, M002, M003, M004
   - Even powers of u,v: M000, M200, M020, M020, M400, M040, M002, M220, M202, M022
   
2. **Sign flip (odd power of original v):**
   - Single u or v: M100→-M010, M010→-M100
   - Odd total (u,v) power with odd v power

