# tenes_simple / tenes_std フェルミオン対応 実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 正方格子・最近接に限定して、`tenes_simple` → `tenes_std` → `tenes` の3段パイプラインで
スピンレスフェルミオンと Hubbard 模型を解けるようにする。

**Architecture:** フェルミオン符号は `tenes_simple` に置く「Fock 構成器」1箇所だけで生成する。
構成器は Jordan–Wigner 規約で2サイト Fock 空間の行列要素を組み、それを rank-4
`op[in1, in2, out1, out2]` に並べ替えて既存スキーマにそのまま載せる。`tenes_std` は
`expm(-τh)` と `parity` の通し役に徹し、2次元ネットワークへの埋め込みで生じる交差符号は
すべて C++ 側の graded 機構が持つ。ツール側の主な追加作業は、この構成器と、
対応範囲外の入力を実際に止めるガード群である。

**Tech Stack:** Python 3(numpy / scipy / toml)、pytest、CMake + CTest、C++17

**Spec:** `docs/superpowers/specs/2026-08-20-fermion-tools-design.md`

## Global Constraints

- 対応範囲は**正方格子・最近接ボンドのみ**。それ以外は明示エラーで止める(黙って通さない)
- 利用者向けのエラー・警告文に `M1` / `M2` を書かない。「現行版では未対応」等にする
- Python は `black`(line-length 88、`tool/pyproject.toml`)、C++ は `clang-format` を通す
- 既存の ctest 28件は全タスクを通じて緑を維持する。ボゾン経路と `test/data/output_*/` の
  golden ファイルには触らない
- 2サイト観測量は fermion では常に `elements` 形式で出力する。`ops = [i, j]` は使わない
- モード順序は `(site1 の spin..., site2 の spin...)`、サイト内は `|up dn> = c†_up c†_dn |0>`
- 局所基底は spinless が `|0>, |1>`(`parity = [0, 1]`)、spinful が
  `|0>, |up>, |dn>, |up dn>`(`parity = [0, 1, 1, 0]`)
- ビルドとテストは `cmake --build out-gcc/build -j 8` と `ctest --test-dir out-gcc/build` を使う
  (ローカルの `gcc` プリセット。macOS なので MPI は OFF)

---

## ファイル構成

| ファイル | 役割 | 変更 |
|---|---|---|
| `tool/tenes_simple.py` | Fock 構成器(モジュール関数)、`Model` の拡張点、フェルミオン2模型、ガード、スキーマ出力 | 変更 |
| `tool/tenes_std.py` | `parity` のパース/検証/出力、fermion 時の入力検証、長距離ボンド拒否 | 変更 |
| `src/iTPS/load_toml.cpp` | `ops` 形式2サイト観測量の拒否 | 変更 |
| `test/python/test_fermion_fock.py` | Fock 構成器の単体テスト(層1・1b) | 新規 |
| `test/python/test_fermion_models.py` | 模型クラス・ガード・スキーマ出力のテスト(層2・3b) | 新規 |
| `test/python/test_tenes_simple.py` | `Model` 拡張点の挙動不変テスト | 追記 |
| `test/python/test_tenes_std.py` | `parity` 通し・長距離拒否のテスト(層3・3b) | 追記 |
| `test/fermion/free_fermion_simple.py.in` | 3段パイプラインの E2E(層4) | 新規 |
| `test/CMakeLists.txt` | E2E テストの登録 | 変更 |
| `sample/08_spinless_fermion_simple/` | simple モードのサンプル | 新規 |
| `docs/sphinx/{ja,en}/file_specification/simple_format.rst` | 模型仕様 | 変更 |
| `NEWS.md` | リリースノート | 変更 |

Fock 構成器を独立モジュールに切らず `tool/tenes_simple.py` のモジュール関数として置くのは、
`tool/` の2スクリプトが `#!<python>` を前置してそのままインストールされる形式であり、
新しい import 先を増やすとインストール手順(`tool/CMakeLists.txt`)まで波及するため。
テストは `sys.path` に `tool/` を挿す既存の流儀(`test/python/test_tenes_simple.py:22-26`)で届く。

---

## Task 1: Fock 構成器

**Files:**
- Modify: `tool/tenes_simple.py`(`dump_op` の直後、94行目付近にモジュール関数を追加)
- Test: `test/python/test_fermion_fock.py`(新規)

**Interfaces:**
- Consumes: なし(numpy のみ)
- Produces:
  - `fermion_modes(nspin: int) -> int`
  - `fock_cop(dagger: bool, mode: int, nmodes: int) -> np.ndarray` — 形状 `(2**nmodes, 2**nmodes)`、軸は `mat[out, in]`
  - `local_index_to_occupation(i: int, nspin: int) -> List[int]` — 長さ `nspin` の占有数リスト
  - `bond_matrix(fock_op: np.ndarray, nspin: int) -> np.ndarray` — 形状 `(d, d, d, d)`、`d = 2**nspin`、軸は `[in1, in2, out1, out2]`
  - `onesite_matrix(fock_op_1site: np.ndarray) -> np.ndarray` — 形状 `(d, d)`、軸は `[in, out]`

設計書 §3 の `onesite_matrix` は `nspin` を第2引数に取っていたが、実装は転置のみで
`nspin` を使わないため落とす(YAGNI)。

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_fermion_fock.py` を新規作成する。冒頭15行の GPL ヘッダは
`test/python/test_tenes_simple.py` の1〜15行目をそのままコピーする。

```python
import os
import sys

import numpy as np
import pytest

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tool")
)

import tenes_simple


def cd(m, M):
    return tenes_simple.fock_cop(True, m, M)


def c(m, M):
    return tenes_simple.fock_cop(False, m, M)


def test_fermion_modes():
    assert tenes_simple.fermion_modes(1) == 2
    assert tenes_simple.fermion_modes(2) == 4


def test_anticommutation_relations():
    M = 4
    dim = 1 << M
    for m in range(M):
        for n in range(M):
            expected = np.eye(dim) if m == n else np.zeros((dim, dim))
            assert np.allclose(cd(m, M) @ c(n, M) + c(n, M) @ cd(m, M), expected)
            assert np.allclose(c(m, M) @ c(n, M) + c(n, M) @ c(m, M), 0.0)
            assert np.allclose(cd(m, M) @ cd(n, M) + cd(n, M) @ cd(m, M), 0.0)


def test_local_index_to_occupation_spinless():
    assert tenes_simple.local_index_to_occupation(0, 1) == [0]
    assert tenes_simple.local_index_to_occupation(1, 1) == [1]


def test_local_index_to_occupation_spinful():
    # |0>, |up>, |dn>, |up dn>  ->  i = n_up + 2 n_dn
    assert tenes_simple.local_index_to_occupation(0, 2) == [0, 0]
    assert tenes_simple.local_index_to_occupation(1, 2) == [1, 0]
    assert tenes_simple.local_index_to_occupation(2, 2) == [0, 1]
    assert tenes_simple.local_index_to_occupation(3, 2) == [1, 1]


def test_bond_matrix_leg_order_is_in_in_out_out():
    # c^dag_{1} c_{2} is not hermitian, so a transposed leg order would show up.
    M = 2
    t = tenes_simple.bond_matrix(cd(0, M) @ c(1, M), 1)
    # <1, 0| c^dag_1 c_2 |0, 1> = +1
    assert t[0, 1, 1, 0] == pytest.approx(1.0)
    # the adjoint element must NOT be here
    assert t[1, 0, 0, 1] == pytest.approx(0.0)


def test_spinless_hopping_matches_the_handwritten_sample():
    # sample/07_spinless_fermion/input.toml, [[observable.twosite]] bond_hamiltonian
    M = 2
    h = -(cd(0, M) @ c(1, M) + cd(1, M) @ c(0, M))
    t = tenes_simple.bond_matrix(h, 1)
    expected = np.zeros((2, 2, 2, 2))
    expected[0, 1, 1, 0] = -1.0
    expected[1, 0, 0, 1] = -1.0
    assert np.allclose(t, expected)


def test_known_sign_for_spinful_hopping():
    # <up dn, 0| c^dag_{1 up} c_{2 up} |dn, up> = -1
    M = 4
    t = tenes_simple.bond_matrix(cd(0, M) @ c(2, M), 2)
    assert t[2, 1, 3, 0] == pytest.approx(-1.0)


def test_doublon_is_cdag_up_cdag_dn_on_the_vacuum():
    # |up dn> = c^dag_up c^dag_dn |0>, so the amplitude is +1 (not -1)
    vac = np.zeros(4)
    vac[0] = 1.0
    state = cd(0, 2) @ (cd(1, 2) @ vac)
    assert state[3] == pytest.approx(1.0)


def test_bond_hamiltonian_is_hermitian_and_parity_conserving():
    M = 4
    parity = [0, 1, 1, 0]
    h = np.zeros((1 << M, 1 << M))
    for s in range(2):
        m1, m2 = s, 2 + s
        h = h - (cd(m1, M) @ c(m2, M) + cd(m2, M) @ c(m1, M))
    assert np.allclose(h, h.conj().T)
    t = tenes_simple.bond_matrix(h, 2)
    for i1, i2, o1, o2 in np.ndindex(t.shape):
        if (parity[i1] ^ parity[i2]) != (parity[o1] ^ parity[o2]):
            assert t[i1, i2, o1, o2] == 0.0


def test_zero_parameters_give_a_zero_bond_matrix():
    M = 2
    h = 0.0 * (cd(0, M) @ c(1, M))
    assert np.allclose(tenes_simple.bond_matrix(h, 1), 0.0)


def test_onesite_matrix_transposes_to_in_out():
    cdag_up = tenes_simple.fock_cop(True, 0, 2)  # mat[out, in]
    m = tenes_simple.onesite_matrix(cdag_up)  # op[in, out]
    assert cdag_up[1, 0] == pytest.approx(1.0)  # <up| c^dag_up |0> = 1
    assert m[0, 1] == pytest.approx(1.0)


def test_onesite_number_operators():
    n_up = tenes_simple.onesite_matrix(cd(0, 2) @ c(0, 2))
    n_dn = tenes_simple.onesite_matrix(cd(1, 2) @ c(1, 2))
    assert np.allclose(np.diag(n_up), [0.0, 1.0, 0.0, 1.0])
    assert np.allclose(np.diag(n_dn), [0.0, 0.0, 1.0, 1.0])
```

- [ ] **Step 2: テストが失敗することを確認する**

```bash
python3 -m pytest test/python/test_fermion_fock.py -v
```

Expected: 全件 FAIL(`AttributeError: module 'tenes_simple' has no attribute 'fermion_modes'`)

- [ ] **Step 3: 最小限の実装を書く**

`tool/tenes_simple.py` の `dump_op` 関数の直後(94行目付近、`Bond = namedtuple(...)` の直前)に
以下を挿入する。

```python
def fermion_modes(nspin: int) -> int:
    """Number of fermion modes in a two-site bond.

    Modes are ordered as (site1 spin...), (site2 spin...).
    """
    return 2 * nspin


def fock_cop(dagger: bool, mode: int, nmodes: int) -> np.ndarray:
    """Creation/annihilation operator including the Jordan-Wigner string.

    The basis is the occupation bit string g, where bit m is mode m and

        |n_0 ... n_{M-1}> = (c^dag_0)^{n_0} ... (c^dag_{M-1})^{n_{M-1}} |0>.

    Moving c^dag_m (or c_m) past the preceding creation operators gives the
    sign (-1)^{sum_{k<m} n_k}.

    The returned matrix follows the NumPy convention ``mat[out, in]``, which is
    the transpose of the ``op[in, out]`` convention used by the one-site
    operators of this module.  ``bond_matrix`` and ``onesite_matrix`` absorb
    the difference.
    """
    dim = 1 << nmodes
    mat = np.zeros((dim, dim))
    for state in range(dim):
        occupied = (state >> mode) & 1
        if dagger == bool(occupied):
            continue
        sign = 1.0
        for k in range(mode):
            if (state >> k) & 1:
                sign = -sign
        mat[state ^ (1 << mode), state] = sign
    return mat


def local_index_to_occupation(i: int, nspin: int) -> List[int]:
    """Local basis index -> occupation numbers (n_up[, n_dn]).

    spinless: i = n.  spinful: i = n_up + 2 n_dn, that is
    |0>, |up>, |dn>, |up dn> for i = 0, 1, 2, 3.  The intra-site order is
    fixed to |up dn> = c^dag_up c^dag_dn |0>.
    """
    return [(i >> s) & 1 for s in range(nspin)]


def bond_matrix(fock_op: np.ndarray, nspin: int) -> np.ndarray:
    """Two-site Fock matrix -> rank-4 op[in1, in2, out1, out2].

    op[i1, i2, o1, o2] = <o1 o2| O |i1 i2> = fock_op[out_global, in_global],
    where the global occupation index is  g = i_site1 + 2**nspin * i_site2
    because the local index is itself the occupation bit pattern.

    The leg order matches what ``tenes_simple`` already emits for a product of
    two one-site operators (``np.einsum("ij,kl -> ikjl", oi, oj)`` with
    ``oi[in, out]``), and matches the plain matrix elements that the C++
    ``wrap_twosite_gate`` / ``wrap_reduced_pair_op`` expect.
    """
    d = 1 << nspin
    op = np.zeros((d, d, d, d), dtype=fock_op.dtype)
    for i1 in range(d):
        for i2 in range(d):
            gin = i1 + d * i2
            for o1 in range(d):
                for o2 in range(d):
                    op[i1, i2, o1, o2] = fock_op[o1 + d * o2, gin]
    return op


def onesite_matrix(fock_op_1site: np.ndarray) -> np.ndarray:
    """One-site Fock matrix mat[out, in] -> op[in, out]."""
    return np.array(fock_op_1site).T
```

`List` は既に `typing` から import 済み(`tool/tenes_simple.py` 冒頭)なので追加 import は不要。

- [ ] **Step 4: テストが通ることを確認する**

```bash
python3 -m pytest test/python/test_fermion_fock.py -v
```

Expected: 全件 PASS

- [ ] **Step 5: 整形してコミットする**

```bash
black --line-length 88 tool/tenes_simple.py test/python/test_fermion_fock.py
git add tool/tenes_simple.py test/python/test_fermion_fock.py
git commit -m "Add the Jordan-Wigner Fock builder to tenes_simple"
```

---

## Task 2: `Model` の拡張点(既存模型は挙動不変)

**Files:**
- Modify: `tool/tenes_simple.py`(`Model` クラス 512-528行、`tenes_simple()` の初期状態出力 1167-1189行、2サイト観測量出力 1267-1295行付近)
- Test: `test/python/test_tenes_simple.py`(追記)

**Interfaces:**
- Consumes: なし
- Produces:
  - `Model.is_fermion: bool` — クラス属性、既定 `False`
  - `Model.parity: List[int]` — インスタンス属性、既定 `[]`(ボゾン模型は空)
  - `Model.twosite_ops_explicit: List[Tuple[str, np.ndarray]]` — インスタンス属性、既定 `[]`。
    各要素は `(name, op)` で `op` は rank-4 `[in1, in2, out1, out2]`
  - `Model.initial_state_vectors(mode: str, num_sublattice: int) -> Optional[np.ndarray]` —
    `None` は「乱数初期化」を意味する。返り値の第0軸は**非 vacancy 副格子を出現順に詰めた添字**

このタスクではフェルミオン模型を一切追加しない。拡張点だけを入れて、既存の出力が
1バイトも変わらないことを確認してから次へ進む。以降の失敗をフェルミオン側に切り分けるため。

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_tenes_simple.py` の末尾に追記する。ファイル冒頭の import に
`import numpy as np` と `import toml` が無ければ追加する。

```python
class TestModelExtensionPoints:
    def test_bosonic_models_are_not_fermionic(self):
        assert tenes_simple.SpinModel({"type": "spin"}).is_fermion is False
        assert tenes_simple.BoseHubbardModel({"type": "boson"}).is_fermion is False

    def test_bosonic_models_have_no_parity_and_no_explicit_twosite_ops(self):
        model = tenes_simple.SpinModel({"type": "spin"})
        assert model.parity == []
        assert model.twosite_ops_explicit == []

    def test_random_mode_returns_none(self):
        model = tenes_simple.SpinModel({"type": "spin"})
        assert model.initial_state_vectors("random", 2) is None

    def test_ferro_mode_repeats_the_first_sublattice(self):
        model = tenes_simple.SpinModel({"type": "spin"})
        st = model.initial_states(2)
        v = model.initial_state_vectors("ferro", 2)
        assert np.allclose(v[0], st[0])
        assert np.allclose(v[1], st[0])

    def test_other_modes_pass_the_pattern_through(self):
        model = tenes_simple.SpinModel({"type": "spin"})
        v = model.initial_state_vectors("antiferro", 2)
        assert np.allclose(v, model.initial_states(2))


def _unitcell_initial_states(std_toml_text):
    parsed = toml.loads(std_toml_text)
    return [u["initial_state"] for u in parsed["tensor"]["unitcell"]]


class TestVacancyInitialStateIndexing:
    """The kagome lattice has a vacancy sublattice; the non-vacancy sublattices
    must keep reading the pattern by their own running index."""

    def _param(self):
        return {
            "parameter": {"general": {}},
            "lattice": {"type": "kagome lattice", "L": 2, "W": 2,
                        "virtual_dim": 2, "initial": "antiferro"},
            "model": {"type": "spin"},
        }

    def test_vacancy_gets_the_scalar_state(self):
        text, _ = tenes_simple.tenes_simple(self._param())
        states = _unitcell_initial_states(text)
        assert states[3] == [1.0]

    def test_nonvacancy_sublattices_follow_the_pattern(self):
        text, lattice = tenes_simple.tenes_simple(self._param())
        states = _unitcell_initial_states(text)
        model = tenes_simple.make_model(self._param())
        pattern = model.initial_states(3)
        for i in range(3):
            assert np.allclose(states[i], pattern[i])
```

- [ ] **Step 2: テストが失敗することを確認する**

```bash
python3 -m pytest test/python/test_tenes_simple.py -v -k "ExtensionPoints or Vacancy"
```

Expected: `TestModelExtensionPoints` は `AttributeError` で FAIL。
`TestVacancyInitialStateIndexing` は現行実装でも PASS する(kagome では vacancy が
最後の副格子なので添字が偶然一致している)。これは**回帰の網**であり、Step 3 の書き換えで
壊れないことを保証する。

- [ ] **Step 3: 拡張点を実装する**

3-a. `Model` クラス(`tool/tenes_simple.py` 512行目付近)の属性宣言と `__init__` を変更する。

```python
class Model(abc.ABC):
    N: int
    onesite_ops: List[np.ndarray]
    onesite_ops_name: List[str]
    twosite_ops: List[Tuple[int, int]]
    twosite_ops_name: List[str]
    twosite_ops_explicit: List[Tuple[str, np.ndarray]]
    params_onesite: Dict[str, Any]  # [neighbor_level][bond_type]
    params_twosite: List[List[Dict[str, Any]]]  # [neighbor_level][bond_type]
    ham_twosites_list: List[List[Tuple[int, int]]]

    is_fermion: bool = False

    def __init__(self):
        self.N = 0
        self.onesite_ops = []
        self.twosite_ops_explicit = []
        self.parity = []
        self.params_onesite = {}
        self.params_twosite = [[]]
        self.ham_twosites_list = [[]]

    def initial_state_vectors(
        self, mode: str, num_sublattice: int
    ) -> Optional[np.ndarray]:
        """Initial product state for each non-vacancy sublattice.

        Returns None when the state should be initialized randomly.  The first
        axis of the returned array is indexed by non-vacancy sublattices in
        order of appearance.
        """
        if mode == "random":
            return None
        st = self.initial_states(num_sublattice)
        if mode == "ferro":
            return np.array([st[0, :] for _ in range(num_sublattice)])
        return st
```

`Optional` が `typing` の import に無ければ追加する。

3-b. `tenes_simple()` の初期状態出力(1167-1189行目付近)を、非 vacancy 用のカウンタを
使う形に置き換える。

```python
    num_sublattice = 0
    for sl in lattice.sublattice:
        if not sl.is_vacancy:
            num_sublattice += 1
    st = model.initial_state_vectors(lattice.initial_states, num_sublattice)
    nonvacancy_index = 0
    for sl in lattice.sublattice:
        ret.append("[[tensor.unitcell]]")
        ret.append("virtual_dim = {}".format(sl.vdim))
        ret.append("index = {}".format(sl.sites))
        if sl.is_vacancy:
            ret.append("physical_dim = {}".format(1))
            ret.append("initial_state = [1.0]")
        else:
            ret.append("physical_dim = {}".format(model.N))
            if st is None:
                state = [0.0]
            else:
                state = st[nonvacancy_index, :]
            nonvacancy_index += 1
            v = ", ".join(map(str, state))
            ret.append("initial_state = [{}]".format(v))
        ret.append("noise = {}".format(lattice.noise))
        ret.append("")
```

3-c. 2サイト観測量の出力ループ(1267-1295行目付近)の直後、`ret.append("")` で終わる
ブロックの後に、明示 rank-4 チャネルの出力を追加する。既存ループが使っている
グループ番号カウンタ `k` をそのまま継続する。

```python
    for name, op in model.twosite_ops_explicit:
        if not (is_complex or np.all(np.isreal(op))):
            continue
        ret.append("[[observable.twosite]]")
        ret.append('name = "{}"'.format(name))
        ret.append("group = {}".format(k))
        k += 1
        ret.append("dim = {}".format([model.N] * 2))
        ret.append('bonds = """')
        for bond in chain(*lattice.bonds[0]):
            ret.append(dumpbond(bond))
        ret.append('"""')
        ret.append('elements = """')
        for line in dump_op(op):
            ret.append(line)
        ret.append('"""')
        ret.append("")
```

- [ ] **Step 4: テストが通ることを確認する**

```bash
python3 -m pytest test/python/test_tenes_simple.py -v
python3 -m pytest test/python -q
```

Expected: 全件 PASS

- [ ] **Step 5: 既存の ctest が緑であることを確認する(このタスクの本題)**

```bash
cmake --build out-gcc/build -j 8
ctest --test-dir out-gcc/build --output-on-failure
```

Expected: 28/28 PASS。1件でも落ちたら**先へ進まず**、拡張点の書き換えを見直す。

- [ ] **Step 6: 整形してコミットする**

```bash
black --line-length 88 tool/tenes_simple.py test/python/test_tenes_simple.py
git add tool/tenes_simple.py test/python/test_tenes_simple.py
git commit -m "Add Model extension points for fermionic models"
```

---

## Task 3: `SpinlessFermionModel` とスキーマ出力

**Files:**
- Modify: `tool/tenes_simple.py`(`BoseHubbardModel` の直後にクラスを追加、`make_model` 1100-1106行、`tenes_simple()` の `[parameter]` 出力 1151-1159行、`[[tensor.unitcell]]` 出力)
- Test: `test/python/test_fermion_models.py`(新規)

**Interfaces:**
- Consumes: Task 1 の `fock_cop` / `bond_matrix` / `onesite_matrix`、Task 2 の `Model.is_fermion` / `Model.parity` / `Model.twosite_ops_explicit` / `Model.initial_state_vectors`
- Produces:
  - `class SpinlessFermionModel(Model)` — `type = "spinless fermion"` で選択される。
    `N = 2`、`parity = [0, 1]`、`is_fermion = True`
  - `tenes_simple()` が `[parameter.general] fermion = true` を自動注入し、
    `[[tensor.unitcell]]` に `parity = [...]` を出力する

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_fermion_models.py` を新規作成する。冒頭15行の GPL ヘッダは
`test/python/test_tenes_simple.py` の1〜15行目をそのままコピーする。

```python
import os
import sys

import numpy as np
import pytest
import toml

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tool")
)

import tenes_simple


def spinless_param(model_extra=None, lattice_extra=None):
    model = {"type": "spinless fermion", "t": 1.0}
    model.update(model_extra or {})
    lattice = {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2}
    lattice.update(lattice_extra or {})
    return {"parameter": {"general": {}}, "lattice": lattice, "model": model}


def std_toml(param):
    text, _ = tenes_simple.tenes_simple(param)
    return toml.loads(text)


class TestSpinlessFermionModel:
    def test_is_selected_by_type(self):
        model = tenes_simple.make_model(spinless_param())
        assert isinstance(model, tenes_simple.SpinlessFermionModel)

    def test_physical_dimension_and_parity(self):
        model = tenes_simple.make_model(spinless_param())
        assert model.N == 2
        assert model.parity == [0, 1]
        assert model.is_fermion is True

    def test_bond_hamiltonian_matches_the_handwritten_sample(self):
        # t = 1, v = 0, mu = 0  ->  sample/07_spinless_fermion/input.toml
        model = tenes_simple.make_model(spinless_param())
        h = model.bondhamiltonian(0, 0, z=4)
        expected = np.zeros((2, 2, 2, 2))
        expected[0, 1, 1, 0] = -1.0
        expected[1, 0, 0, 1] = -1.0
        assert np.allclose(h, expected)

    def test_chemical_potential_is_split_over_the_bonds(self):
        model = tenes_simple.make_model(spinless_param({"mu": 2.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        # -mu/z * (n1 + n2), z = 4  ->  -0.5 on each occupied site
        assert h[1, 0, 1, 0] == pytest.approx(-0.5)
        assert h[0, 1, 0, 1] == pytest.approx(-0.5)
        assert h[1, 1, 1, 1] == pytest.approx(-1.0)

    def test_nearest_neighbour_repulsion(self):
        model = tenes_simple.make_model(spinless_param({"v": 3.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        assert h[1, 1, 1, 1] == pytest.approx(3.0)

    def test_onesite_observable_is_the_density(self):
        model = tenes_simple.make_model(spinless_param())
        assert model.onesite_ops_name == ["n"]
        assert np.allclose(model.onesite_ops[0], np.diag([0.0, 1.0]))

    def test_hopping_is_an_explicit_rank4_observable(self):
        model = tenes_simple.make_model(spinless_param())
        names = [name for name, _ in model.twosite_ops_explicit]
        assert "hopping" in names
        op = dict(model.twosite_ops_explicit)["hopping"]
        assert op.shape == (2, 2, 2, 2)
        assert op[0, 1, 1, 0] == pytest.approx(1.0)
        assert op[1, 0, 0, 1] == pytest.approx(1.0)


class TestSpinlessFermionSchema:
    def test_fermion_flag_is_injected(self):
        assert std_toml(spinless_param())["parameter"]["general"]["fermion"] is True

    def test_explicit_fermion_false_is_rejected(self):
        param = spinless_param()
        param["parameter"]["general"]["fermion"] = False
        with pytest.raises(RuntimeError, match="fermion"):
            tenes_simple.tenes_simple(param)

    def test_parity_is_emitted_for_every_unitcell(self):
        parsed = std_toml(spinless_param())
        for ucell in parsed["tensor"]["unitcell"]:
            assert ucell["parity"] == [0, 1]

    def test_bosonic_models_do_not_emit_parity_or_fermion(self):
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "spin"},
        }
        parsed = std_toml(param)
        assert "fermion" not in parsed["parameter"]["general"]
        for ucell in parsed["tensor"]["unitcell"]:
            assert "parity" not in ucell

    def test_no_twosite_observable_uses_the_ops_form(self):
        text, _ = tenes_simple.tenes_simple(spinless_param())
        assert "ops = " not in text

    def test_vacuum_initial_state(self):
        parsed = std_toml(spinless_param(lattice_extra={"initial": "vacuum"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert np.allclose(ucell["initial_state"], [1.0, 0.0])

    def test_random_initial_state_stays_scalar(self):
        parsed = std_toml(spinless_param(lattice_extra={"initial": "random"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert ucell["initial_state"] == [0.0]

    @pytest.mark.parametrize("mode", ["ferro", "antiferro", "full", "cdw"])
    def test_unsupported_initial_states_are_rejected(self, mode):
        with pytest.raises(RuntimeError):
            tenes_simple.tenes_simple(spinless_param(lattice_extra={"initial": mode}))
```

- [ ] **Step 2: テストが失敗することを確認する**

```bash
python3 -m pytest test/python/test_fermion_models.py -v
```

Expected: 全件 FAIL(`RuntimeError: Unknown model type: spinless fermion`)

- [ ] **Step 3: モデルクラスとスキーマ出力を実装する**

3-a. `BoseHubbardModel` の直後(980行目付近、`def make_lattice` の手前)に追加する。

```python
class SpinlessFermionModel(Model):
    """Spinless fermions on the square lattice (nearest neighbour only).

    Local basis |0>, |1> with fermion parity [0, 1].

        H_bond = -t (c^dag_1 c_2 + h.c.) + V n_1 n_2 - (mu/z) (n_1 + n_2)
    """

    is_fermion = True

    def __init__(self, param: Dict[str, Any]):
        super().__init__()
        self.N = 2
        self.nspin = 1
        self.parity = [0, 1]

        nmodes = 1
        n_op = fock_cop(True, 0, nmodes) @ fock_cop(False, 0, nmodes)
        self.onesite_ops = [onesite_matrix(n_op)]
        self.onesite_ops_name = ["n"]

        self.twosite_ops = [(0, 0)]
        self.twosite_ops_name = ["nn"]

        M = fermion_modes(self.nspin)
        hop = fock_cop(True, 0, M) @ fock_cop(False, 1, M)
        hop = hop + fock_cop(True, 1, M) @ fock_cop(False, 0, M)
        self.twosite_ops_explicit = [("hopping", bond_matrix(hop, self.nspin))]

        self.read_params(param)

    def read_params(self, modelparam: Dict[str, Any]) -> None:
        ret_onesite: Dict[str, Any] = {}
        ret_twosite: List[List[Dict[str, Any]]] = [
            [{}, {}, {}],  # 1st neighbors
            [{}, {}, {}],  # 2nd neighbors
            [{}, {}, {}],  # 3rd neighbors
        ]

        repat = re.compile("^([tv])([012]?)('{0,2})$")
        for key in modelparam.keys():
            if key in ("type", "mu"):
                continue
            ma = repat.match(key)
            if not ma:
                msg = "Unknown keyname {}".format(key)
                raise RuntimeError(msg)
            gr = ma.groups()
            types = [int(gr[1])] if gr[1] else [0, 1, 2]
            n = len(gr[2])
            for typ in types:
                if gr[0] in ret_twosite[n][typ]:
                    raise RuntimeError("{} is defined twice".format(key))
                ret_twosite[n][typ][gr[0]] = modelparam[key]

        ret_onesite["mu"] = modelparam.get("mu", 0.0)
        self.params_onesite = ret_onesite
        self.params_twosite = ret_twosite
        # ham_twosites_list is rebuilt by Model.sort_ham_groups(), which
        # hamiltonians() calls before using it; do not set it here.

    def initial_states(self, num_sublattice: int) -> np.ndarray:
        ret = np.zeros((num_sublattice, self.N))
        ret[:, 0] = 1.0
        return ret

    def initial_state_vectors(
        self, mode: str, num_sublattice: int
    ) -> Optional[np.ndarray]:
        if mode == "random":
            return None
        if mode == "vacuum":
            return self.initial_states(num_sublattice)
        msg = 'initial = "{}" is not available for spinless fermions'.format(mode)
        msg += '; use "random" or "vacuum".'
        msg += " A product state with an odd-parity site (such as |1>) cannot be"
        msg += " built, because TeNeS puts the state vector on virtual index 0"
        msg += " (even) and the total leg parity of the site tensor would be odd."
        raise RuntimeError(msg)

    def model_sitehamiltonian(self, params_onesite: Dict) -> np.ndarray:
        return np.zeros((self.N, self.N))

    def model_bondhamiltonian(
        self,
        z: int,
        use_onesite_hamiltonian: bool,
        params_onesite: Dict,
        params_twosite: Dict,
    ) -> np.ndarray:
        t = params_twosite.get("t", 0.0)
        V = params_twosite.get("v", 0.0)
        mu = params_onesite.get("mu", 0.0)

        M = fermion_modes(self.nspin)

        def cd(m):
            return fock_cop(True, m, M)

        def cop(m):
            return fock_cop(False, m, M)

        n1 = cd(0) @ cop(0)
        n2 = cd(1) @ cop(1)
        h = -t * (cd(0) @ cop(1) + cd(1) @ cop(0))
        h = h + V * (n1 @ n2)
        h = h - (mu / z) * (n1 + n2)
        return bond_matrix(h, self.nspin)
```

`re` は既に import 済み。

3-b. `make_model`(1100-1106行目)に分岐を足す。`"square lattice"` と同じ流儀で
`startswith` を使う。

```python
    modelparam = param["model"]
    model: Model
    if modelparam["type"] == "spin":
        model = SpinModel(modelparam)
    elif modelparam["type"] == "boson":
        model = BoseHubbardModel(modelparam)
    elif modelparam["type"].startswith("spinless"):
        model = SpinlessFermionModel(modelparam)
    else:
        msg = "Unknown model type: {}".format(modelparam["type"])
        raise RuntimeError(msg)
    return model
```

3-c. `tenes_simple()` の `[parameter]` 出力(1151-1159行目)の直前に `fermion` の
自動注入を入れる。

```python
    ret = []
    ret.append("[parameter]")
    pparam = param["parameter"]
    if model.is_fermion:
        general = pparam.setdefault("general", {})
        if general.get("fermion", True) is False:
            msg = "parameter.general.fermion = false conflicts with the model type"
            msg += ' "{}", which is fermionic.'.format(param["model"]["type"])
            raise RuntimeError(msg)
        general["fermion"] = True
    for name in ("general", "simple_update", "full_update", "ctm", "random"):
        ...
```

3-d. `[[tensor.unitcell]]` の出力(Task 2 で書き換えたブロック)に `parity` を足す。
vacancy 側は `parity = [0]`。

```python
        if sl.is_vacancy:
            ret.append("physical_dim = {}".format(1))
            if model.is_fermion:
                ret.append("parity = [0]")
            ret.append("initial_state = [1.0]")
        else:
            ret.append("physical_dim = {}".format(model.N))
            if model.is_fermion:
                ret.append("parity = {}".format(model.parity))
            ...
```

- [ ] **Step 4: テストが通ることを確認する**

```bash
python3 -m pytest test/python/test_fermion_models.py -v
python3 -m pytest test/python -q
```

Expected: 全件 PASS

- [ ] **Step 5: 整形してコミットする**

```bash
black --line-length 88 tool/tenes_simple.py test/python/test_fermion_models.py
git add tool/tenes_simple.py test/python/test_fermion_models.py
git commit -m "Add the spinless fermion model to tenes_simple"
```

---

## Task 4: 対応範囲外の入力を止めるガード

**Files:**
- Modify: `tool/tenes_simple.py`(`tenes_simple()` の `model = make_model(param)` 直後、1148行目付近)
- Test: `test/python/test_fermion_models.py`(追記)

**Interfaces:**
- Consumes: Task 3 の `SpinlessFermionModel`、`Model.is_fermion`
- Produces: `_check_fermion_scope(param, lattice, model) -> None` — モジュール内のヘルパー。
  範囲外なら `RuntimeError`

設計書 §6.3。方針を書くだけでは止まらないので、ここで実際の拒否を入れる。

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_fermion_models.py` の末尾に追記する。

```python
class TestFermionScopeGuards:
    @pytest.mark.parametrize(
        "latname", ["honeycomb lattice", "triangular lattice", "kagome lattice"]
    )
    def test_non_square_lattices_are_rejected(self, latname):
        param = spinless_param(lattice_extra={"type": latname})
        with pytest.raises(RuntimeError, match="square"):
            tenes_simple.tenes_simple(param)

    def test_square_lattice_is_accepted(self):
        tenes_simple.tenes_simple(spinless_param())

    # In read_params the digit is the BOND TYPE and the number of primes is
    # the NEIGHBOUR LEVEL, so t1 / t2 are still nearest neighbour and only the
    # primed keys go beyond it.
    @pytest.mark.parametrize("key", ["t'", "t''", "v'", "v''"])
    def test_beyond_nearest_neighbour_parameters_are_rejected(self, key):
        param = spinless_param({key: 0.5})
        with pytest.raises(RuntimeError, match="nearest"):
            tenes_simple.tenes_simple(param)

    def test_zero_valued_far_neighbour_parameters_are_accepted(self):
        tenes_simple.tenes_simple(spinless_param({"t'": 0.0}))

    def test_bond_type_variants_of_the_first_neighbour_are_accepted(self):
        tenes_simple.tenes_simple(spinless_param({"t0": 1.0}))

    def test_correlation_is_rejected(self):
        param = spinless_param()
        param["correlation"] = {"r_max": 5, "operators": [[0, 0]]}
        with pytest.raises(RuntimeError, match="correlation"):
            tenes_simple.tenes_simple(param)

    def test_correlation_length_is_rejected(self):
        param = spinless_param()
        param["correlation_length"] = {"measure": True}
        with pytest.raises(RuntimeError, match="correlation_length"):
            tenes_simple.tenes_simple(param)

    def test_bosonic_models_are_untouched_by_the_guards(self):
        param = {
            "parameter": {"general": {}},
            "lattice": {"type": "kagome lattice", "L": 2, "W": 2, "virtual_dim": 2},
            "model": {"type": "spin", "j": 1.0, "j'": 0.5},
            "correlation": {"r_max": 3, "operators": [[0, 0]]},
        }
        tenes_simple.tenes_simple(param)

    def test_the_message_does_not_mention_the_internal_milestone(self):
        param = spinless_param(lattice_extra={"type": "honeycomb lattice"})
        with pytest.raises(RuntimeError) as excinfo:
            tenes_simple.tenes_simple(param)
        assert "M1" not in str(excinfo.value)
```

- [ ] **Step 2: テストが失敗することを確認する**

```bash
python3 -m pytest test/python/test_fermion_models.py -v -k ScopeGuards
```

Expected: 拒否系のテストが FAIL(`RuntimeError` が上がらない)

- [ ] **Step 3: ガードを実装する**

3-a. `tenes_simple()` の `model = make_model(param)`(1148行目)と
`hams = hamiltonians(...)`(1149行目)の間に呼び出しを挿入する。

```python
    param = lower_dict(param)
    lattice = make_lattice(param)
    model = make_model(param)
    _check_fermion_scope(param, lattice, model)
    hams = hamiltonians(lattice, model, use_onesite_hamiltonian)
```

3-b. `def tenes_simple(` の直前にヘルパーを定義する。

```python
def _check_fermion_scope(
    param: MutableMapping[str, Any], lattice: Lattice, model: Model
) -> None:
    """Reject inputs outside the supported fermionic scope.

    The current version supports the square lattice with nearest-neighbour
    bonds only.  Nothing else is silently converted; every unsupported input
    stops here with a reason.
    """
    if not model.is_fermion:
        return

    scope = (
        "the fermion support in this version covers the square lattice with"
        " nearest-neighbour bonds only"
    )

    if not isinstance(lattice, SquareLattice):
        msg = 'lattice type "{}" is not available for fermionic models; {}.'.format(
            param["lattice"]["type"], scope
        )
        raise RuntimeError(msg)

    for n, per_level in enumerate(model.params_twosite):
        if n == 0:
            continue
        for typ, params in enumerate(per_level):
            for name, value in params.items():
                if value != 0.0:
                    msg = "{} = {} is a {}-neighbour term; {}.".format(
                        name, value, n + 1, scope
                    )
                    raise RuntimeError(msg)

    if "correlation" in param:
        msg = "[correlation] is not available for fermionic models in this version"
        msg += "; remove the section."
        raise RuntimeError(msg)

    if "correlation_length" in param:
        msg = "[correlation_length] is not available for fermionic models in this"
        msg += " version; remove the section. The transfer-matrix correlation"
        msg += " length is not fermion-aware, and the solver would silently"
        msg += " disable it."
        raise RuntimeError(msg)
```

- [ ] **Step 4: テストが通ることを確認する**

```bash
python3 -m pytest test/python/test_fermion_models.py -v
python3 -m pytest test/python -q
```

Expected: 全件 PASS

- [ ] **Step 5: 整形してコミットする**

```bash
black --line-length 88 tool/tenes_simple.py test/python/test_fermion_models.py
git add tool/tenes_simple.py test/python/test_fermion_models.py
git commit -m "Reject inputs outside the supported fermionic scope"
```

---

## Task 5: `tenes_std` の `parity` 通しと長距離ボンド拒否

**Files:**
- Modify: `tool/tenes_std.py`(`LocalTensor.__init__` 275-290行、`LocalTensor.check` 292-312行、`tenes_std.Model.__init__` の 1177行目直前、`to_toml` の unitcell 出力 1207-1215行)
- Test: `test/python/test_tenes_std.py`(追記)

**Interfaces:**
- Consumes: Task 3 が出力する `std.toml`(`fermion = true`、`parity`)
- Produces: `input.toml` に `parity` が現れる。fermion 時の検証が
  時間発展演算子の生成**前**に走る

**最重要**: `make_evolution_twosite()` は `nhops > 1` のとき、中間サイトに恒等演算子を
テンソル積してから SVD で最近接ゲート列へ分解する。この分解は JW string を置かないので
フェルミオンでは誤りであり、しかも分解後の各ゲートは最近接なので C++ 側の距離ガードには
到達しない。**黙って間違った数値が出る唯一の経路**なので、ここで止める。

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_tenes_std.py` の末尾に追記する。ファイル冒頭の import に
`import numpy as np`、`import pytest`、`import toml`、`import io` が無ければ追加する。

```python
def _fermion_std_param(bond="0 1 0", extra_general=None):
    general = {"fermion": True}
    general.update(extra_general or {})
    return {
        "parameter": {
            "general": general,
            "simple_update": {"tau": 0.01, "num_step": 10},
            "ctm": {"dimension": 4},
        },
        "tensor": {
            "L_sub": [2, 2],
            "unitcell": [
                {
                    "index": [],
                    "physical_dim": 2,
                    "virtual_dim": [2, 2, 2, 2],
                    "parity": [0, 1],
                    "initial_state": [0.0],
                    "noise": 1.0,
                }
            ],
        },
        "hamiltonian": [
            {
                "dim": [2, 2],
                "bonds": bond,
                "elements": "0 1 1 0 -1.0 0.0\n1 0 0 1 -1.0 0.0\n",
            }
        ],
        "observable": {"onesite": [], "twosite": []},
    }


def _to_toml_text(inp):
    buf = io.StringIO()
    inp.to_toml(buf)
    return buf.getvalue()


class TestFermionParityPassthrough:
    def test_parity_is_parsed(self):
        inp = tenes_std.Model(_fermion_std_param())
        assert inp.unitcell.sites[0].parity == [0, 1]

    def test_parity_is_written_out(self):
        inp = tenes_std.Model(_fermion_std_param())
        parsed = toml.loads(_to_toml_text(inp))
        assert parsed["tensor"]["unitcell"][0]["parity"] == [0, 1]

    def test_bosonic_input_has_no_parity(self):
        param = _fermion_std_param()
        del param["parameter"]["general"]["fermion"]
        del param["tensor"]["unitcell"][0]["parity"]
        inp = tenes_std.Model(param)
        parsed = toml.loads(_to_toml_text(inp))
        assert "parity" not in parsed["tensor"]["unitcell"][0]

    def test_missing_parity_is_rejected_in_fermion_mode(self):
        param = _fermion_std_param()
        del param["tensor"]["unitcell"][0]["parity"]
        with pytest.raises(RuntimeError, match="parity"):
            tenes_std.Model(param)

    def test_parity_with_wrong_length_is_rejected(self):
        param = _fermion_std_param()
        param["tensor"]["unitcell"][0]["parity"] = [0, 1, 0]
        with pytest.raises(RuntimeError, match="parity"):
            tenes_std.Model(param)

    def test_parity_with_invalid_values_is_rejected(self):
        param = _fermion_std_param()
        param["tensor"]["unitcell"][0]["parity"] = [0, 2]
        with pytest.raises(RuntimeError, match="parity"):
            tenes_std.Model(param)


class TestFermionLongDistanceBond:
    def test_a_distance_two_bond_is_rejected(self):
        # (dx, dy) = (2, 0) needs two hops, which would be split into a chain
        # of nearest-neighbour gates without a Jordan-Wigner string.
        param = _fermion_std_param(bond="0 2 0")
        with pytest.raises(RuntimeError, match="nearest"):
            tenes_std.Model(param)

    def test_a_diagonal_bond_is_rejected(self):
        param = _fermion_std_param(bond="0 1 1")
        with pytest.raises(RuntimeError, match="nearest"):
            tenes_std.Model(param)

    def test_a_nearest_neighbour_bond_is_accepted(self):
        tenes_std.Model(_fermion_std_param(bond="0 1 0"))

    def test_bosonic_input_still_allows_long_bonds(self):
        param = _fermion_std_param(bond="0 2 0")
        del param["parameter"]["general"]["fermion"]
        del param["tensor"]["unitcell"][0]["parity"]
        tenes_std.Model(param)


class TestFermionObservableForm:
    def test_ops_form_twosite_observable_is_rejected(self):
        param = _fermion_std_param()
        param["observable"]["onesite"] = [
            {"group": 0, "name": "n", "sites": [], "dim": 2, "elements": "1 1 1.0 0.0\n"}
        ]
        param["observable"]["twosite"] = [
            {"group": 1, "name": "nn", "bonds": "0 1 0", "dim": [2, 2], "ops": [0, 0]}
        ]
        with pytest.raises(RuntimeError, match="elements"):
            tenes_std.Model(param)
```

- [ ] **Step 2: テストが失敗することを確認する**

```bash
python3 -m pytest test/python/test_tenes_std.py -v -k "Fermion"
```

Expected: `AttributeError: 'LocalTensor' object has no attribute 'parity'` などで FAIL

- [ ] **Step 3: `parity` の保持・検証・出力と長距離ボンド拒否を実装する**

3-a. `LocalTensor` のクラス属性宣言(`tool/tenes_std.py` 261-264行目付近)に追加する。

```python
    phys_dim: int
    virtual_dim: List[int]
    parity: Optional[List[int]]
```

3-b. `LocalTensor.__init__`(284行目付近)でパースする。

```python
        self.phys_dim = tensor_dict["physical_dim"]

        virtual_dim = tensor_dict["virtual_dim"]
        if isinstance(virtual_dim, int):
            self.virtual_dim = [virtual_dim] * 4
        else:
            self.virtual_dim = virtual_dim
        self.parity = tensor_dict.get("parity", None)
        self.check()
```

3-c. `LocalTensor.check`(292行目付近)の末尾に検証を足す。

```python
        if self.parity is not None:
            if not (
                isinstance(self.parity, list)
                and len(self.parity) == self.phys_dim
                and all(p in (0, 1) for p in self.parity)
            ):
                msg = "parity must be a list of {} values, each 0 or 1".format(
                    self.phys_dim
                )
                raise RuntimeError(msg)
```

3-d. `to_toml` の unitcell 出力(1207-1215行目付近)に足す。

```python
            f.write("physical_dim = {}\n".format(ucell["physical_dim"]))
            f.write("virtual_dim = {}\n".format(ucell["virtual_dim"]))
            if "parity" in ucell:
                f.write("parity = {}\n".format(ucell["parity"]))
```

3-e. `tenes_std.Model.__init__` の `self.simple_updates = []`(1177行目)の**直前**に
fermion 検証を挿入する。時間発展演算子を作る前に止める必要がある。

```python
        self._check_fermion_input()

        self.simple_updates = []
        self.full_updates = []
```

3-f. `tenes_std.Model` のメソッドとして追加する(`to_toml` の直前)。

```python
    def _check_fermion_input(self) -> None:
        """Reject inputs that the fermionic path cannot handle.

        Runs before the evolution operators are built, because the long-bond
        decomposition below would otherwise produce nearest-neighbour gates
        that look valid to the solver but carry no Jordan-Wigner string.
        """
        general = self.parameter.get("general", {})
        if not general.get("fermion", False):
            return

        for i, site in enumerate(self.unitcell.sites):
            if site.parity is None:
                msg = "tensor.unitcell[{}] has no parity, which is required".format(i)
                msg += " in fermion mode."
                raise RuntimeError(msg)

        for obs in self.twobodies:
            if getattr(obs, "ops", None):
                msg = 'two-site observable "{}" uses the ops form, which is not'.format(
                    obs.name
                )
                msg += " available in fermion mode; write it out as elements."
                raise RuntimeError(msg)

        for ham in self.hamiltonians:
            if isinstance(ham, SiteOperator):
                continue
            nhops = len(self.graph.make_path(ham.bond))
            if nhops != 1:
                source = ham.bond.source_site
                msg = "a bond term connecting site {} and (dx, dy) = ({}, {})".format(
                    source, ham.bond.dx, ham.bond.dy
                )
                msg += " is not a nearest-neighbour bond; the fermion support in"
                msg += " this version covers nearest-neighbour bonds only."
                msg += " Longer bonds would be split into a chain of gates"
                msg += " without the Jordan-Wigner string."
                raise RuntimeError(msg)
```

`self.twobodies` の各要素が `ops` 属性を持つかは実装依存なので `getattr` で読む。
`Optional` が `typing` の import に無ければ追加する。

- [ ] **Step 4: テストが通ることを確認する**

```bash
python3 -m pytest test/python/test_tenes_std.py -v
python3 -m pytest test/python -q
```

Expected: 全件 PASS

- [ ] **Step 5: 3段パイプラインが `input.toml` まで通ることを手で確認する**

リポジトリのルートを `REPO` に入れて実行する。

```bash
REPO=$(git rev-parse --show-toplevel)
WORK=$(mktemp -d)
cd "$WORK"
cat > simple.toml <<'EOT'
[parameter.general]
is_real = true
[parameter.simple_update]
tau = 0.01
num_step = 10
[parameter.ctm]
dimension = 8
[lattice]
type = "square lattice"
L = 2
W = 2
virtual_dim = 2
initial = "random"
noise = 1.0
[model]
type = "spinless fermion"
t = 1.0
EOT
python3 "$REPO/tool/tenes_simple.py" simple.toml -o std.toml
python3 "$REPO/tool/tenes_std.py" std.toml -o input.toml
grep -n "fermion\|parity\|ops = " input.toml
```

Expected: `fermion = True`(または `true`)と `parity = [0, 1]` が現れ、`ops = ` は現れない。

- [ ] **Step 6: 整形してコミットする**

```bash
black --line-length 88 tool/tenes_std.py test/python/test_tenes_std.py
git add tool/tenes_std.py test/python/test_tenes_std.py
git commit -m "Carry parity through tenes_std and reject long fermionic bonds"
```

---

## Task 6: τ 展開テスト(`tenes_simple` と `tenes_std` を結ぶ)

**Files:**
- Test: `test/python/test_fermion_models.py`(追記)

**Interfaces:**
- Consumes: Task 3 の `SpinlessFermionModel`、Task 5 の `tenes_std.tenes_std.Model`
- Produces: なし(テストのみ)

`tenes_simple` が作った `h` と `tenes_std` が作ったゲートを結ぶ唯一の検査。
τ → 0 で `(expm(-τh) - I)/(-τ) → h` に戻ることを確認する。

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_fermion_models.py` の末尾に追記する。冒頭の import に
`import tenes_std` を追加する(`sys.path` は既に `tool/` を含む)。

```python
def gate_and_hamiltonian(param, tau):
    """Run the tenes_simple -> tenes_std pipeline and return (gate, h_bond).

    Every nearest-neighbour bond of the uniform square lattice carries the same
    gate, so simple_updates[0] is representative.
    """
    param = dict(param)
    param["parameter"] = {
        "general": {"is_real": True},
        "simple_update": {"tau": tau, "num_step": 1},
        "ctm": {"dimension": 4},
    }
    text, _ = tenes_simple.tenes_simple(param)
    inp = tenes_std.Model(toml.loads(text))
    gate = inp.simple_updates[0].elements
    ham = None
    for h in inp.hamiltonians:
        elements = getattr(h, "elements", None)
        if elements is not None and elements.ndim == 4:
            ham = elements
            break
    assert ham is not None, "no two-site hamiltonian was produced"
    return gate, ham


def assert_gate_expands_to_hamiltonian(param, tau=1.0e-6):
    gate, ham = gate_and_hamiltonian(param, tau)
    d = ham.shape[0]
    # identity on [in1, in2, out1, out2]
    identity = np.einsum("ik,jl->ijkl", np.eye(d), np.eye(d))
    recovered = (gate - identity) / (-tau)
    assert np.allclose(recovered, ham, atol=1e-4)


def test_spinless_gate_expands_to_the_hamiltonian():
    assert_gate_expands_to_hamiltonian(spinless_param({"t": 1.0, "v": 0.7, "mu": 0.3}))
```

- [ ] **Step 2: テストを走らせる**

```bash
python3 -m pytest test/python/test_fermion_models.py -v -k gate_expands
```

Expected: PASS。FAIL する場合、脚順か符号がどこかでずれている。
`gate` と `ham` の shape、`inp.simple_updates[0]` の中身を印字して原因を切り分ける。
このテストは**新しい実装を要求しない**——Task 1〜5 が正しければ通る。通らなければ
そこに実バグがある。

- [ ] **Step 3: コミットする**

```bash
black --line-length 88 test/python/test_fermion_models.py
git add test/python/test_fermion_models.py
git commit -m "Pin that the fermionic gate expands back to the bond Hamiltonian"
```

---

## Task 7: E2E テスト(3段パイプライン)

**Files:**
- Create: `test/fermion/free_fermion_simple.py.in`
- Modify: `test/CMakeLists.txt`(100-107行目付近、`FreeFermion` の登録の直後)

**Interfaces:**
- Consumes: Task 3〜5 の全経路
- Produces: ctest エントリ `FreeFermionSimple`

**配位数について**: `SquareLattice.zs[0] = [2, 2, 0]` なので、最近接ボンドを1つの
ハミルトニアン群にまとめたときの `ztotal` は 2 + 2 + 0 = 4 になる。これが
`model_bondhamiltonian` に渡る `z` であり、正方格子の配位数と一致する。
テストで `bondhamiltonian(0, 0, z=4)` と書いているのはこのため。

パラメータは手書きサンプル `sample/07_spinless_fermion/input.toml` と揃える
(t = 1、mu = 0、D = 2、χ = 8、tau = 0.01、num_step = 1000、seed = 11、noise = 1.0)。
サンプルの実測は E = −0.7328、n = 0.5000、実行時間 45 秒。

- [ ] **Step 1: テストスクリプトを書く**

`test/fermion/free_fermion_simple.py.in` を新規作成する。GPL ヘッダは
`test/fermion/free_fermion.py.in` の冒頭からコピーし、CMake の置換変数
(`@TENES_PYTHON_EXECUTABLE@` 等)の使い方も同ファイルに倣う。

```python
from __future__ import print_function

import os
import subprocess
import sys
from os.path import join

import toml

PYTHON = "@TENES_PYTHON_EXECUTABLE@"
TOOL_DIR = "@CMAKE_SOURCE_DIR@/tool"
TENES = "@CMAKE_BINARY_DIR@/src/tenes"
WORK_DIR = "@CMAKE_CURRENT_BINARY_DIR@/fermion/free_fermion_simple"

SIMPLE_TOML = """
[parameter.general]
is_real = true
output = "output"

[parameter.simple_update]
tau = 0.01
num_step = 1000

[parameter.ctm]
dimension = 8
convergence_epsilon = 1.0e-8
iteration_max = 100

[parameter.random]
seed = 11

[lattice]
type = "square lattice"
L = 2
W = 2
virtual_dim = 2
initial = "random"
noise = 1.0

[model]
type = "spinless fermion"
t = 1.0
mu = 0.0
"""


def run(cmd):
    print("run:", " ".join(cmd), flush=True)
    subprocess.check_call(cmd, cwd=WORK_DIR)


def run_tenes(inputfile):
    # same launcher construction as test/fermion/free_fermion.py.in
    cmd = []
    if "@MPIEXEC@":
        cmd.append("@MPIEXEC@")
        cmd.extend("@MPIEXEC_PREFLAGS@".split())
        cmd.append("@MPIEXEC_NUMPROC_FLAG@")
        cmd.append("1")
        cmd.extend("@MPIEXEC_POSTFLAGS@".split())
    cmd.append(TENES)
    cmd.append(inputfile)
    run(cmd)


def read_density(outdir):
    # density.dat lines look like
    #   Energy           = -7.32844895452057554e-01  0.00000000000000000e+00
    # so the real part is the FIRST word after "=", not the last.
    values = {}
    with open(join(outdir, "density.dat")) as f:
        for line in f:
            if "=" not in line:
                continue
            name, raw = line.split("=", 1)
            words = raw.split()
            values[name.strip()] = float(words[0])
    return values


def main():
    if not os.path.isdir(WORK_DIR):
        os.makedirs(WORK_DIR)
    with open(join(WORK_DIR, "simple.toml"), "w") as f:
        f.write(SIMPLE_TOML)

    run([PYTHON, join(TOOL_DIR, "tenes_simple.py"), "simple.toml", "-o", "std.toml"])
    run([PYTHON, join(TOOL_DIR, "tenes_std.py"), "std.toml", "-o", "input.toml"])

    generated = toml.load(join(WORK_DIR, "input.toml"))
    assert generated["parameter"]["general"]["fermion"] is True, "fermion flag missing"
    for ucell in generated["tensor"]["unitcell"]:
        assert ucell["parity"] == [0, 1], "parity missing or wrong"
    with open(join(WORK_DIR, "input.toml")) as f:
        assert "ops = " not in f.read(), "two-site observable in ops form"

    run_tenes("input.toml")

    values = read_density(join(WORK_DIR, "output"))
    energy = values["Energy"]
    density = values["n"]

    print("Energy =", energy, " n =", density, flush=True)
    assert abs(density - 0.5) < 5.0e-3, \
        "density {} is not at half filling".format(density)
    # Exact free-fermion ground state energy per site: -0.81056.
    # D = 2 is a coarse approximation; the handwritten expert-mode sample
    # sample/07_spinless_fermion gives -7.32844895452057554e-01 with the same
    # t, mu, D, chi, tau, num_step, seed and noise.
    assert -0.85 < energy < -0.70, \
        "energy {} is out of the expected range".format(energy)
    print("OK", flush=True)


main()
```

- [ ] **Step 2: CMake に登録する**

`test/CMakeLists.txt` の `FreeFermion` 登録ブロック(100-107行目付近)の直後に追加する。

```cmake
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fermion/free_fermion_simple.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_simple.py @ONLY)
add_test(NAME FreeFermionSimple
         COMMAND ${TENES_PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/fermion/free_fermion_simple.py)
set_tests_properties(FreeFermionSimple PROPERTIES TIMEOUT 600)
```

- [ ] **Step 3: 走らせて実測値を得る**

```bash
cmake --build out-gcc/build -j 8
ctest --test-dir out-gcc/build -R FreeFermionSimple --output-on-failure
```

Expected: PASS。標準出力に `Energy = ... n = ...` が出る。実行時間は1〜2分以内。

- [ ] **Step 4: 実測値で回帰を締める**

手書きサンプル `sample/07_spinless_fermion/output/density.dat` の値は
`Energy = -7.32844895452057554e-01`、`n = 4.99999684464887739e-01` である。
パイプラインが同じ入力を作れていれば、Step 3 の実測値はこれに一致するはずである。

一致していたら、緩い範囲判定を回帰ピンに置き換える。

```python
    # matches the handwritten expert-mode sample bit-for-bit as far as the
    # printed precision goes; see sample/07_spinless_fermion/output/density.dat
    E_REF = -7.32844895452057554e-01
    assert abs(energy - E_REF) < 1.0e-2 * abs(E_REF), \
        "energy {} moved away from the reference {}".format(energy, E_REF)
```

一致していなければ、パイプラインがサンプルと違う入力を作っている。
`diff <(sort input.toml) <(sort ../../../sample/07_spinless_fermion/input.toml)` の要領で
生成された `input.toml` をサンプルと突き合わせ、原因を突き止めてから先へ進む。
差が evolution ゲートの並び順だけで数値が一致しているなら、実測値を `E_REF` にして
その旨をコメントに書く。

- [ ] **Step 5: 再度走らせて通ることを確認し、コミットする**

```bash
ctest --test-dir out-gcc/build -R FreeFermionSimple --output-on-failure
git add test/fermion/free_fermion_simple.py.in test/CMakeLists.txt
git commit -m "Add an end-to-end test for the fermionic three-stage pipeline"
```

---

## Task 8: `HubbardModel`

**Files:**
- Modify: `tool/tenes_simple.py`(`SpinlessFermionModel` の直後にクラスを追加、`make_model` に分岐)
- Test: `test/python/test_fermion_models.py`(追記)

**Interfaces:**
- Consumes: Task 1 の Fock 構成器、Task 2〜4 の拡張点とガード
- Produces: `class HubbardModel(Model)` — `type = "hubbard"` で選択される。
  `N = 4`、`parity = [0, 1, 1, 0]`、`is_fermion = True`

d = 2 では reduced-pair blob の「素ロード」と「入力+出力 swap」が数保存演算子で縮退するため、
Task 3 の spinless 一致テストでは d = 4 の規約を捕まえられない。ここで独立に固定する。

- [ ] **Step 1: 失敗するテストを書く**

`test/python/test_fermion_models.py` の末尾に追記する。

```python
def hubbard_param(model_extra=None, lattice_extra=None):
    model = {"type": "hubbard", "t": 1.0}
    model.update(model_extra or {})
    lattice = {"type": "square lattice", "L": 2, "W": 2, "virtual_dim": 2}
    lattice.update(lattice_extra or {})
    return {"parameter": {"general": {}}, "lattice": lattice, "model": model}


class TestHubbardModel:
    def test_is_selected_by_type(self):
        model = tenes_simple.make_model(hubbard_param())
        assert isinstance(model, tenes_simple.HubbardModel)

    def test_physical_dimension_and_parity(self):
        model = tenes_simple.make_model(hubbard_param())
        assert model.N == 4
        assert model.parity == [0, 1, 1, 0]
        assert model.is_fermion is True

    def test_hopping_carries_the_jordan_wigner_sign(self):
        # <up dn, 0| H |dn, up> comes from -t c^dag_{1 up} c_{2 up}, whose
        # matrix element is -1, so the entry is +t.
        model = tenes_simple.make_model(hubbard_param())
        h = model.bondhamiltonian(0, 0, z=4)
        assert h[2, 1, 3, 0] == pytest.approx(1.0)

    def test_hubbard_u_appears_on_doubly_occupied_sites(self):
        model = tenes_simple.make_model(hubbard_param({"u": 8.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        # U/z on each site, z = 4  ->  site1 doubly occupied, site2 empty
        assert h[3, 0, 3, 0] == pytest.approx(2.0)
        assert h[3, 3, 3, 3] == pytest.approx(4.0)

    def test_bond_hamiltonian_is_parity_conserving(self):
        model = tenes_simple.make_model(hubbard_param({"u": 4.0, "v": 1.0, "mu": 2.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        parity = model.parity
        for i1, i2, o1, o2 in np.ndindex(h.shape):
            if (parity[i1] ^ parity[i2]) != (parity[o1] ^ parity[o2]):
                assert h[i1, i2, o1, o2] == 0.0

    def test_bond_hamiltonian_is_hermitian(self):
        model = tenes_simple.make_model(hubbard_param({"u": 4.0, "v": 1.0, "mu": 2.0}))
        h = model.bondhamiltonian(0, 0, z=4)
        # rows are (in1, in2), columns are (out1, out2)
        m = h.reshape(16, 16)
        assert np.allclose(m, m.conj().T)

    def test_onesite_observables(self):
        model = tenes_simple.make_model(hubbard_param())
        assert model.onesite_ops_name == [
            "n", "n_up", "n_dn", "Sz", "doublon", "holon"
        ]
        ops = dict(zip(model.onesite_ops_name, model.onesite_ops))
        assert np.allclose(np.diag(ops["n"]), [0.0, 1.0, 1.0, 2.0])
        assert np.allclose(np.diag(ops["n_up"]), [0.0, 1.0, 0.0, 1.0])
        assert np.allclose(np.diag(ops["n_dn"]), [0.0, 0.0, 1.0, 1.0])
        assert np.allclose(np.diag(ops["Sz"]), [0.0, 0.5, -0.5, 0.0])
        assert np.allclose(np.diag(ops["doublon"]), [0.0, 0.0, 0.0, 1.0])
        assert np.allclose(np.diag(ops["holon"]), [1.0, 0.0, 0.0, 0.0])

    def test_every_onesite_observable_is_parity_even(self):
        model = tenes_simple.make_model(hubbard_param())
        parity = model.parity
        for op in model.onesite_ops:
            for i, o in np.ndindex(op.shape):
                if parity[i] != parity[o]:
                    assert op[i, o] == 0.0


class TestHubbardSchema:
    def test_parity_is_emitted(self):
        parsed = std_toml(hubbard_param())
        for ucell in parsed["tensor"]["unitcell"]:
            assert ucell["parity"] == [0, 1, 1, 0]

    def test_vacuum_initial_state(self):
        parsed = std_toml(hubbard_param(lattice_extra={"initial": "vacuum"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert np.allclose(ucell["initial_state"], [1.0, 0.0, 0.0, 0.0])

    def test_full_initial_state(self):
        parsed = std_toml(hubbard_param(lattice_extra={"initial": "full"}))
        for ucell in parsed["tensor"]["unitcell"]:
            assert np.allclose(ucell["initial_state"], [0.0, 0.0, 0.0, 1.0])

    def test_cdw_initial_state_alternates(self):
        parsed = std_toml(hubbard_param(lattice_extra={"initial": "cdw"}))
        states = [u["initial_state"] for u in parsed["tensor"]["unitcell"]]
        assert np.allclose(states[0], [1.0, 0.0, 0.0, 0.0])
        assert np.allclose(states[1], [0.0, 0.0, 0.0, 1.0])

    @pytest.mark.parametrize("mode", ["ferro", "antiferro", "neel"])
    def test_unsupported_initial_states_are_rejected(self, mode):
        with pytest.raises(RuntimeError):
            tenes_simple.tenes_simple(hubbard_param(lattice_extra={"initial": mode}))

    def test_guards_apply_to_hubbard_too(self):
        with pytest.raises(RuntimeError, match="square"):
            tenes_simple.tenes_simple(
                hubbard_param(lattice_extra={"type": "triangular lattice"})
            )


def test_hubbard_gate_expands_to_the_hamiltonian():
    assert_gate_expands_to_hamiltonian(
        hubbard_param({"t": 1.0, "u": 4.0, "v": 0.5, "mu": 1.0})
    )
```

- [ ] **Step 2: テストが失敗することを確認する**

```bash
python3 -m pytest test/python/test_fermion_models.py -v -k "Hubbard or hubbard"
```

Expected: `RuntimeError: Unknown model type: hubbard` で FAIL

- [ ] **Step 3: `HubbardModel` を実装する**

`SpinlessFermionModel` の直後に追加する。

```python
class HubbardModel(Model):
    """Fermionic Hubbard model on the square lattice (nearest neighbour only).

    Local basis |0>, |up>, |dn>, |up dn> with fermion parity [0, 1, 1, 0], and
    the intra-site order fixed to |up dn> = c^dag_up c^dag_dn |0>.

        H_bond = -t sum_s (c^dag_{1s} c_{2s} + h.c.) + V n_1 n_2
                 + (1/z) [ U (n_up n_dn)_1 + U (n_up n_dn)_2
                           - mu (n_1 + n_2) - h (Sz_1 + Sz_2) ]

    Note this is the *fermionic* Hubbard model.  The Bose-Hubbard model is
    ``type = "boson"``.
    """

    is_fermion = True

    def __init__(self, param: Dict[str, Any]):
        super().__init__()
        self.N = 4
        self.nspin = 2
        self.parity = [0, 1, 1, 0]

        ns = self.nspin

        def cd1(m):
            return fock_cop(True, m, ns)

        def c1(m):
            return fock_cop(False, m, ns)

        n_up = cd1(0) @ c1(0)
        n_dn = cd1(1) @ c1(1)
        n_tot = n_up + n_dn
        sz = 0.5 * (n_up - n_dn)
        doublon = n_up @ n_dn
        holon = (np.eye(4) - n_up) @ (np.eye(4) - n_dn)

        self.onesite_ops = [
            onesite_matrix(n_tot),
            onesite_matrix(n_up),
            onesite_matrix(n_dn),
            onesite_matrix(sz),
            onesite_matrix(doublon),
            onesite_matrix(holon),
        ]
        self.onesite_ops_name = ["n", "n_up", "n_dn", "Sz", "doublon", "holon"]

        # (index into onesite_ops) pairs: nn = n*n, SzSz = Sz*Sz
        self.twosite_ops = [(0, 0), (3, 3)]
        self.twosite_ops_name = ["nn", "SzSz"]

        M = fermion_modes(self.nspin)
        hop = np.zeros((1 << M, 1 << M))
        for s in range(self.nspin):
            m1, m2 = s, self.nspin + s
            hop = hop + fock_cop(True, m1, M) @ fock_cop(False, m2, M)
            hop = hop + fock_cop(True, m2, M) @ fock_cop(False, m1, M)
        self.twosite_ops_explicit = [("hopping", bond_matrix(hop, self.nspin))]

        self.read_params(param)

    def read_params(self, modelparam: Dict[str, Any]) -> None:
        ret_onesite: Dict[str, Any] = {}
        ret_twosite: List[List[Dict[str, Any]]] = [
            [{}, {}, {}],  # 1st neighbors
            [{}, {}, {}],  # 2nd neighbors
            [{}, {}, {}],  # 3rd neighbors
        ]

        repat = re.compile("^([tv])([012]?)('{0,2})$")
        for key in modelparam.keys():
            if key in ("type", "mu", "u", "h"):
                continue
            ma = repat.match(key)
            if not ma:
                raise RuntimeError("Unknown keyname {}".format(key))
            gr = ma.groups()
            types = [int(gr[1])] if gr[1] else [0, 1, 2]
            n = len(gr[2])
            for typ in types:
                if gr[0] in ret_twosite[n][typ]:
                    raise RuntimeError("{} is defined twice".format(key))
                ret_twosite[n][typ][gr[0]] = modelparam[key]

        ret_onesite["mu"] = modelparam.get("mu", 0.0)
        ret_onesite["u"] = modelparam.get("u", 0.0)
        ret_onesite["h"] = modelparam.get("h", 0.0)
        self.params_onesite = ret_onesite
        self.params_twosite = ret_twosite
        # see the note in SpinlessFermionModel.read_params

    def initial_states(self, num_sublattice: int) -> np.ndarray:
        ret = np.zeros((num_sublattice, self.N))
        ret[:, 0] = 1.0
        return ret

    def initial_state_vectors(
        self, mode: str, num_sublattice: int
    ) -> Optional[np.ndarray]:
        if mode == "random":
            return None
        ret = np.zeros((num_sublattice, self.N))
        if mode == "vacuum":
            ret[:, 0] = 1.0
            return ret
        if mode == "full":
            ret[:, 3] = 1.0
            return ret
        if mode == "cdw":
            for i in range(num_sublattice):
                ret[i, 3 if i % 2 else 0] = 1.0
            return ret
        msg = 'initial = "{}" is not available for the hubbard model'.format(mode)
        msg += '; use "random", "vacuum", "full", or "cdw".'
        msg += " A product state with an odd-parity site (such as |up>) cannot"
        msg += " be built, because TeNeS puts the state vector on virtual index"
        msg += " 0 (even) and the total leg parity of the site tensor would be"
        msg += " odd."
        raise RuntimeError(msg)

    def model_sitehamiltonian(self, params_onesite: Dict) -> np.ndarray:
        return np.zeros((self.N, self.N))

    def model_bondhamiltonian(
        self,
        z: int,
        use_onesite_hamiltonian: bool,
        params_onesite: Dict,
        params_twosite: Dict,
    ) -> np.ndarray:
        t = params_twosite.get("t", 0.0)
        V = params_twosite.get("v", 0.0)
        mu = params_onesite.get("mu", 0.0)
        U = params_onesite.get("u", 0.0)
        hz = params_onesite.get("h", 0.0)

        ns = self.nspin
        M = fermion_modes(ns)

        def cd(m):
            return fock_cop(True, m, M)

        def cop(m):
            return fock_cop(False, m, M)

        def number(site, spin):
            m = site * ns + spin
            return cd(m) @ cop(m)

        h = np.zeros((1 << M, 1 << M))
        for s in range(ns):
            m1, m2 = s, ns + s
            h = h - t * (cd(m1) @ cop(m2) + cd(m2) @ cop(m1))

        n1 = number(0, 0) + number(0, 1)
        n2 = number(1, 0) + number(1, 1)
        h = h + V * (n1 @ n2)

        for site in (0, 1):
            nu, nd = number(site, 0), number(site, 1)
            h = h + (U / z) * (nu @ nd)
            h = h - (mu / z) * (nu + nd)
            h = h - (hz / z) * 0.5 * (nu - nd)

        return bond_matrix(h, ns)
```

`make_model` に分岐を足す。

```python
    elif modelparam["type"].startswith("spinless"):
        model = SpinlessFermionModel(modelparam)
    elif modelparam["type"] == "hubbard":
        model = HubbardModel(modelparam)
```

- [ ] **Step 4: テストが通ることを確認する**

```bash
python3 -m pytest test/python/test_fermion_models.py -v
python3 -m pytest test/python -q
```

Expected: 全件 PASS

- [ ] **Step 5: 整形してコミットする**

```bash
black --line-length 88 tool/tenes_simple.py test/python/test_fermion_models.py
git add tool/tenes_simple.py test/python/test_fermion_models.py
git commit -m "Add the fermionic Hubbard model to tenes_simple"
```

---

## Task 9: C++ 側の `ops` 形式2サイト観測量の拒否

**Files:**
- Modify: `src/iTPS/load_toml.cpp`(668-682行目、`validate_fermion_constraints` の
  2サイト観測量ループ)
- Test: `test/test_input.cpp` は触らない。手で `input.toml` を作って確認する

**Interfaces:**
- Consumes: なし
- Produces: fermion モードで `ops = [i, j]` 形式の2サイト観測量が `input_error` になる

ツール側の拒否(Task 5)だけでは `input.toml` を直書きする利用者を止められないため、
独立の必須修正である。

- [ ] **Step 1: 現状を再現する入力を作る**

```bash
REPO=$(git rev-parse --show-toplevel)
WORK=$(mktemp -d)
cd "$WORK"
cp "$REPO/sample/07_spinless_fermion/input.toml" .
```

`input.toml` の `[[observable.twosite]]` ブロック(57-74行目)を、`elements` の代わりに
`ops` を使う形に書き換える。

```toml
[[observable.twosite]]
name = "nn"
group = 1
bonds = """
0 1 0
"""
dim = [2, 2]
ops = [0, 0]
```

- [ ] **Step 2: 現状の挙動を確認する**

```bash
"$REPO/out-gcc/build/src/tenes" input.toml
```

Expected: 入力エラーにならず、実行に入る(そのまま完走するか実行時エラーになる)。
どちらにせよ「意図の伝わる入力エラー」にはなっていない。

- [ ] **Step 3: ガードを追加する**

`src/iTPS/load_toml.cpp` の2サイト観測量ループ(668-682行目)を書き換える。

```cpp
  for (const auto &op : twosite_operators) {
    if (!op.dx.empty() &&
        !is_nearest_neighbor_displacement(op.dx[0], op.dy[0])) {
      throw_fermion_guard("distance-2-or-longer two-site operators");
    }
    if (!op.ops_indices.empty()) {
      throw_fermion_guard(
          "two-site observables in the ops form; write them out as elements");
    }
    const int site1 = lattice.other(op.source_site, op.dx[0], op.dy[0]);
    if (has_odd_tensor_element(op.op,
                               two_site_parity(peps_parameters.phys_parity,
                                               op.source_site, site1))) {
      throw_fermion_guard("parity-odd two-site operators");
    }
  }
```

`ops_indices` が非空のとき `op.op` は空なので、パリティ検査より**前**に拒否する必要がある。

- [ ] **Step 4: ガードが発火することを確認する**

```bash
cmake --build "$REPO/out-gcc/build" -j 8
cd "$WORK"
"$REPO/out-gcc/build/src/tenes" input.toml
```

Expected: `fermion mode in this version does not support two-site observables in
the ops form; write them out as elements; disable fermion mode or remove this
setting` というメッセージで停止する。

- [ ] **Step 5: 元のサンプルが通ることと ctest 全件を確認する**

```bash
cd "$REPO/sample/07_spinless_fermion" && "$REPO/out-gcc/build/src/tenes" input.toml
cd "$REPO" && ctest --test-dir out-gcc/build --output-on-failure
```

Expected: サンプルは完走、ctest は全件 PASS

- [ ] **Step 6: 整形してコミットする**

```bash
clang-format -i src/iTPS/load_toml.cpp
git add src/iTPS/load_toml.cpp
git commit -m "Reject two-site observables in ops form in fermion mode"
```

---

## Task 10: ドキュメント・サンプル・NEWS

**Files:**
- Create: `sample/08_spinless_fermion_simple/simple.toml`、`sample/08_spinless_fermion_simple/README.md`
- Modify: `docs/sphinx/ja/file_specification/simple_format.rst`、`docs/sphinx/en/file_specification/simple_format.rst`、`NEWS.md`

**Interfaces:**
- Consumes: Task 3〜8 の全機能
- Produces: なし

- [ ] **Step 1: simple モードのサンプルを作る**

`sample/08_spinless_fermion_simple/simple.toml` に Task 7 の `SIMPLE_TOML` と同じ内容を書く。
`sample/08_spinless_fermion_simple/README.md` には以下を書く。

- 何を計算するか(t = 1、μ = 0 のスピンレス自由フェルミオン、正方格子、半充填)
- 実行方法の3行(`tenes_simple simple.toml` → `tenes_std std.toml` → `tenes input.toml`)
- 期待される結果(D = 2、χ = 8 で E ≈ −0.7328、n = 0.5。厳密値は −0.81056)
- 既存の `07_spinless_fermion` は expert モード(`input.toml` 直書き)の例として残っており、
  同じ物理系を扱っていること
- 現行版の制約(正方格子・最近接のみ、相関関数・相関長・full update・有限温度は未対応)

- [ ] **Step 2: サンプルが動くことを確認する**

```bash
cd sample/08_spinless_fermion_simple
python3 ../../tool/tenes_simple.py simple.toml -o std.toml
python3 ../../tool/tenes_std.py std.toml -o input.toml
../../out-gcc/build/src/tenes input.toml
```

Expected: Task 7 と同じ値。確認後、生成物(`std.toml`、`input.toml`、`output/`)は
コミットしない(`.gitignore` を確認し、必要なら削除する)。

- [ ] **Step 3: Sphinx ドキュメントに模型を追加する**

`docs/sphinx/ja/file_specification/simple_format.rst` の `[model]` 節に、
既存の `"spin"` / `"boson"` と同じ体裁で2つ追加する。

- `type = "spinless fermion"`: パラメータ表(`t`、`v`、`mu`)、局所基底 `|0>, |1>`、
  1サイト観測量 `n`、2サイト観測量 `hopping`、`nn`
- `type = "hubbard"`: パラメータ表(`t`、`u`、`v`、`mu`、`h`)、局所基底
  `|0>, |up>, |dn>, |up dn>`、サイト内順序 `|up dn> = c†_up c†_dn |0>`、
  1サイト観測量 `n`、`n_up`、`n_dn`、`Sz`、`doublon`、`holon`、
  2サイト観測量 `hopping`、`nn`、`SzSz`。
  **これはフェルミオン Hubbard 模型であり、Bose-Hubbard 模型は `"boson"` である**旨を明記

`[lattice] initial` の説明に、フェルミオン模型で取りうる値
(`"random"`、`"vacuum"`、Hubbard のみ `"full"`、`"cdw"`)と、
`"ferro"` / `"antiferro"` が使えない理由(パリティ奇のプロダクト状態は作れない)を追記する。

制約(正方格子・最近接のみ)も明記する。

`docs/sphinx/en/file_specification/simple_format.rst` に同じ内容の英語版を書く。

- [ ] **Step 4: ドキュメントがビルドできることを確認する**

```bash
cmake --preset docs && cmake --build --preset docs
```

Expected: 警告なしでビルドが通る。プリセットが無ければ
`sphinx-build -b html docs/sphinx/ja out-docs/ja` で代替する。

- [ ] **Step 5: NEWS.md に追記する**

未リリース節に以下を1項目として書く。

```markdown
- `tenes_simple` now supports two fermionic models on the square lattice with
  nearest-neighbour bonds: `type = "spinless fermion"` and `type = "hubbard"`
  (the fermionic Hubbard model; the Bose-Hubbard model remains `type = "boson"`).
  `tenes_std` carries the `parity` metadata through to `input.toml`.
```

- [ ] **Step 6: コミットする**

```bash
git add sample/08_spinless_fermion_simple docs/sphinx NEWS.md
git commit -m "Document the fermionic models in tenes_simple"
```

- [ ] **Step 7: 最終確認**

```bash
cmake --build out-gcc/build -j 8
ctest --test-dir out-gcc/build --output-on-failure
python3 -m pytest test/python -q
black --check --line-length 88 tool/ test/python/
git status
```

Expected: ctest 29/29 PASS(既存28件 + `FreeFermionSimple`)、pytest 全件 PASS、
`black --check` が通り、作業ツリーがクリーン。
