#!/bin/bash
# Golden-file self-tests for sgclinalg. Values pinned against GAP 4.15.1:
#   A := [[1,1,0,0],[0,1,1,0],[1,0,1,0]]*Z(2);
#   SolutionMat(A, [1,0,1,0]*Z(2)) = (1,1,0)          -> "sol 110"
#   NullspaceMat(A)                = [(1,1,1)]
#   BaseMat(A) spans {(1100),(0110)}; our canonical RREF form is
#   (1,0,1,0),(0,1,1,0)   (relaxed parity: span-equal, deterministic)
set -e
cd "$(dirname "$0")"
BIN=$(cd ../.. && pwd)/bin/sgclinalg
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() { echo "FAIL: $1"; exit 1; }

cat > A.sms <<'EOF'
3 4 2
1 1 1
1 2 1
2 2 1
2 3 1
3 1 1
3 3 1
0 0 0
EOF

$BIN rank2 A.sms out.txt
grep -qx 'rank 2' out.txt || fail "rank2"

$BIN nullsp2 A.sms out.sms
grep -qx '1 3 2' <(head -1 out.sms) || fail "nullsp2 header"
grep -q '^1 1 1$' out.sms && grep -q '^1 2 1$' out.sms && grep -q '^1 3 1$' out.sms || fail "nullsp2 vector"

cat > B.sms <<'EOF'
2 4 2
1 1 1
1 3 1
2 4 1
0 0 0
EOF
$BIN solve2 A.sms B.sms out.txt
[ "$(sed -n 1p out.txt)" = "sol 110" ] || fail "solve2 row1"
[ "$(sed -n 2p out.txt)" = "fail" ]    || fail "solve2 row2 (e4 not in rowspace)"

$BIN rowbasis2 A.sms out.sms
# canonical RREF: rows (1,0,1,0) and (0,1,1,0)
grep -qx '2 4 2' <(head -1 out.sms) || fail "rowbasis2 header"
grep -q '^1 1 1$' out.sms && grep -q '^1 3 1$' out.sms || fail "rowbasis2 row1"
grep -q '^2 2 1$' out.sms && grep -q '^2 3 1$' out.sms || fail "rowbasis2 row2"
grep -q '^1 2 1$' out.sms && fail "rowbasis2 not fully reduced" || true

# Entries with even accumulated value must vanish mod 2 (duplicate triplets sum).
cat > DUP.sms <<'EOF'
1 2 2
1 1 1
1 1 1
1 2 1
0 0 0
EOF
$BIN rank2 DUP.sms out.txt
grep -qx 'rank 1' out.txt || fail "duplicate-triplet mod-2 accumulation"

# Zero matrix edge case: rank 0, nullspace = identity.
cat > Z.sms <<'EOF'
2 3 2
0 0 0
EOF
$BIN rank2 Z.sms out.txt
grep -qx 'rank 0' out.txt || fail "zero matrix rank"
$BIN nullsp2 Z.sms out.sms
grep -qx '2 2 2' <(head -1 out.sms) || fail "zero matrix nullspace dim"

# --- Task 2 ops ---

$BIN semiech2 A.sms out.txt
grep -qx 'nvec 2 ncols 4 nrowsA 3' <(head -1 out.txt) || fail "semiech2 header"
[ "$(grep -c '^vec '   out.txt)" = 2 ] || fail "semiech2 vec count"
[ "$(grep -c '^coeff ' out.txt)" = 2 ] || fail "semiech2 coeff count"
grep -q '^heads 1 2 0 0$' out.txt || fail "semiech2 heads"

# steinitz2: complement of span{(1100)} inside span{(1100),(0110)}.
# GAP BaseSteinitzVectors gives factorspace [(0,1,1,0)]; our sifting agrees here.
cat > K.sms <<'EOF'
2 4 2
1 1 1
1 2 1
2 2 1
2 3 1
0 0 0
EOF
cat > I.sms <<'EOF'
1 4 2
1 1 1
1 2 1
0 0 0
EOF
$BIN steinitz2 K.sms I.sms out.sms
grep -qx '1 4 2' <(head -1 out.sms) || fail "steinitz2 header"
grep -q '^1 2 1$' out.sms && grep -q '^1 3 1$' out.sms || fail "steinitz2 vector"

# steinitz2 with empty I: complement is a basis of rowspace(K) itself.
cat > I0.sms <<'EOF'
0 4 2
0 0 0
EOF
$BIN steinitz2 K.sms I0.sms out.sms
grep -qx '2 4 2' <(head -1 out.sms) || fail "steinitz2 empty-I"

echo ALL OK
