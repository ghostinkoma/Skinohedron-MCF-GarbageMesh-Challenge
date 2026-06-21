# 00a · Prologue — How the Model Was Searched For / 数理モデルへの模索の足跡

*A prologue, not a result. This chapter records **why** the project could take the
steps it did — the branch points, the choices not taken, and where each choice was
later settled — before the mathematics is stated. It introduces no new claim of its
own; it points to the chapter and the representative formula where each judgement
was verified. In the spirit of 空論を重ねない, even the story of the search is tied
back to something checked.*

*結果ではなく、序章である。本章は数理モデルを述べる前に、なぜこのプロジェクトthatこの
歩みを取れたか——分岐点、選ばなかった道、そして各選択thatのちにどこで決着したか——を
記録する。新たな主張は持たない。各判断that検証された章と代表式へ案内するだけである。
空論を重ねないという原則に従い、模索の物語さえ、確かめられた何かに繋ぎ留める。*

---

## 1. The starting question / 出発点の問い

**EN.** The project did not begin with an operator. It began with a stubborn
question: *can a simulation survive on a mesh that is not ideal?* The working name
— a "garbage-mesh challenge" on a skin-like polyhedron — was literal. Most clean
derivations assume a nice mesh; the wager here was to insist on physics that stays
correct even when the mesh is merely *equally spaced*, not optimal. That insistence
forced every later choice to be defended by a check, because on an imperfect mesh
intuition is not enough.

**JP.** このプロジェクトは作用素から始まったのではない。執拗な問いから始まった——
*理想的でないメッシュの上で、シミュレーションは成立するか*。仮称「ガベージメッシュ・
チャレンジ」は文字どおりの意味だった。多くの綺麗な導出は良質なメッシュを前提とするthat、
ここでの賭けは、メッシュthat最適でなく単に*等間隔*であっても正しくあり続ける物理を要求
することだった。この要求thatのちのすべての選択を「検証で defend する」ことへ追い込んだ。
不完全なメッシュの上では、直感だけでは足りないからである。

---

## 2. First instinct, and the negative result that redirected it / 最初の直感と、それを
転回させた負の結果

**EN.** The first instinct was to read wave propagation as *scattering* on a
simplicial complex — an S-parameter view, where energy hops between nodes
(`01`–`01d`). It was natural and it ran. But a quiet test exposed a defect: the
naive uniform-node wave operator had a **direction-dependent speed**. Measuring the
anisotropy
$$\Delta c/c = \frac{\max_i c(\theta_i) - \min_i c(\theta_i)}{\mathrm{mean}_i\, c(\theta_i)}$$
showed it did not vanish the way a faithful operator's should. This was the first
fork, and the decision was to **reject** the scattered-node operator rather than
paper over it — the same discipline that later kept negative results in the record.

*Settled in:* `01e` (geometric consistency), formula $\Delta c/c$.

**JP.** 最初の直感は、波の伝播を単体複体上の*散乱*として読むことだった——エネルギーthat
ノード間を飛び移る S パラメータ的な見方（`01`–`01d`）。自然で、実際に動いた。だthat静かな
テストthat欠陥を露わにした。素朴な一様ノードの波作用素は、**速度that方向に依存**していた。
異方性
$$\Delta c/c = \frac{\max_i c(\theta_i) - \min_i c(\theta_i)}{\mathrm{mean}_i\, c(\theta_i)}$$
を測ると、忠実な作用素なら消えるべきものthat消えなかった。これthat最初の分岐で、判断は
散乱ノード作用素を糊塗せず**退ける**ことだった——のちに負の結果を記録に残し続けた、その
同じ流儀である。

*決着:* `01e`（幾何的整合性）、式 $\Delta c/c$。

---

## 3. The pivot: conduction, done the right way / 転回点——「正しいやり方」の伝導

**EN.** Rejecting the scattered node forced a better question: what *is* the
conductance between two nodes? The tempting answer — "proportional to contact area"
— is right only when it is the **cotangent / finite-element conductance**, not a
naive graph weight. Re-derived that way, the operator passed the tests the first
one failed: the lowest non-trivial eigenvalue hit the analytic
$$\lambda_1 = \pi^2 = 9.8696,$$
and a two-material interface reproduced the exact contact temperature
$$T = \frac{k_2}{k_1 + k_2}.$$
This was the moment the foundation became trustworthy: the heat route, not the
scatter route, was the integrating one.

*Settled in:* `02` (heat-conduction route), formulas $\lambda_1=\pi^2$,
$T=k_2/(k_1+k_2)$.

**JP.** 散乱ノードを退けたことthat、より良い問いを強いた。二ノード間の伝導とは*何*か。
魅力的な答え——「接触面積に比例」——thatが正しいのは、それthat**コタンジェント／有限要素
の伝導**であるときだけで、素朴なグラフ重みではない。そう導き直すと、作用素は最初のものthat
落ちたテストを通った。最小の非自明固有値that解析値
$$\lambda_1 = \pi^2 = 9.8696$$
に当たり、二材界面thatが厳密な接触温度
$$T = \frac{k_2}{k_1 + k_2}$$
を再現した。ここで土台thatが信頼に足るものになった。統合する経路は散乱ではなく、熱で
あった。

*決着:* `02`（熱伝導の経路）、式 $\lambda_1=\pi^2$, $T=k_2/(k_1+k_2)$。

---

## 4. The bet on one operator / 一つの作用素への賭け

**EN.** With a trustworthy conductance in hand, a larger pattern came into view.
The heat operator and the wave operator were not two things — they were one,
$$L = M^{-1}K,$$
read at two time-structures: heat as decay along its spectrum, waves as oscillation
on the *same* spectrum,
$$\frac{\partial T}{\partial t} = -L\,T \;\;\leftrightarrow\;\; M\ddot p = -K p.$$
This was the central bet of the whole notebook: that *one* geometric operator could
carry every scalar physics that followed. It was a bet precisely because it could
have failed — and the rule was that it would only be believed once each physics was
checked against an exact answer.

*Settled in:* `03` (one operator, two physics), formula $L=M^{-1}K$.

**JP.** 信頼できる伝導を手にすると、より大きなパターンthat見えてきた。熱の作用素と波の
作用素は二つのものではなく、一つ——
$$L = M^{-1}K$$
——を二つの時間構造で読んだものだった。熱はそのスペクトルに沿う減衰、波は*同じ*スペクトル
上の振動として、
$$\frac{\partial T}{\partial t} = -L\,T \;\;\leftrightarrow\;\; M\ddot p = -K p.$$
これthatノート全体の中心的な賭けだった——*一つ*の幾何作用素that、続くすべてのスカラー物理を
担えるという賭け。失敗しうるからこそ賭けであり、規則は「各物理that厳密解と照合されて
初めて信じる」ことだった。

*決着:* `03`（一つの作用素、二つの物理）、式 $L=M^{-1}K$。

---

## 5. What the search method became / 模索that方法になった

**EN.** By this point the *way of searching* had crystallised into the project's
standing method, applied to everything afterward (dynamics and boundaries in `03b`,
pressure in `04`, mesh transform in `05`/`5b`):

1. state a hypothesis as a formula;
2. find the **exact solution** it must reproduce;
3. write code that asserts the match to machine precision;
4. advance only on PASS — and keep the FAILs (the scattered node, the sphere
   obstruction, the un-removable sign-flips) as part of the record.

The prologue's real subject is this method, born from the very first negative
result. The mathematics that follows is downstream of it.

**JP.** この段階で、*模索のやり方*そのものthatプロジェクトの恒常的な方法へ結晶していた。
以後のすべて（`03b` の力学と境界、`04` の圧力、`05`/`5b` のメッシュ変形）に適用される。

1. 仮説を式として述べる；
2. それthat再現すべき**厳密解**を見つける；
3. 機械精度での一致を assert するコードを書く；
4. PASS でのみ前進する——そして FAIL（散乱ノード、球の障害、消せない符号反転）を記録の
一部として残す。

序章の真の主題はこの方法であり、それは最初の負の結果から生まれた。続く数学は、その下流に
ある。

---

## 6. Why step into V3 now / なぜ今 V3 へ踏み出すか

**EN.** The model has now shown two things that, together, justify the next step.
It **unifies** (one `L` for heat, wave, and pressure — `03`, `04`), and it
**transforms** (the same physics carried onto torus and sphere — `05`), with the
honest cost of geometric deformation measured rather than hidden (`5b`). A framework
that both unifies and transforms is exactly one that earns a harder question. V3 is
that harder question — the move beyond static and linear scalar fields toward the
next layer of physics. This prologue does not yet state V3's model; it records that
the ground under it was laid deliberately, checked at every step, so that the step
into V3 is taken from solid footing rather than hope.

*The trajectory in one line:* a stubborn question about bad meshes → a rejected
first instinct → conduction done correctly → one operator → a verified method →
unification and transform → the readiness to ask more.

**JP.** モデルは今、次の一歩を正当化する二つのことを示した。**統合する**（熱・波・圧力に
一つの `L`——`03`, `04`）こと、そして**変形する**（同じ物理thatトーラスと球へ運ばれる——
`05`）こと。幾何変形の正直な代償は、隠さず測られている（`5b`）。統合し、かつ変形する枠組み
こそ、より難しい問いに値する。V3 はその難しい問いである——静的・線形なスカラー場を超え、
次の物理の層へ向かう一歩。本序章はまだ V3 のモデルを述べない。述べるのは、その足場that
意図的に、各段で確かめられて築かれたということ——だから V3 への一歩は、希望からではなく、
固い足場から踏み出される。

*一行での軌跡:* 悪いメッシュへの執拗な問い → 退けた最初の直感 → 正しく行った伝導 →
一つの作用素 → 検証された方法 → 統合と変形 → さらに問う準備。

---

*Continues into V3. The verified chapters this prologue points to — `01e`, `02`,
`03`, `03b`, `04`, `05`/`5b` — and the synthesis `00` are where the claims live;
here only the footsteps are recorded.*

*V3 へ続く。本序章that指す検証済みの章——`01e`, `02`, `03`, `03b`, `04`, `05`/`5b`——
と総括 `00` に主張は宿る。ここに記したのは足跡のみである。*
